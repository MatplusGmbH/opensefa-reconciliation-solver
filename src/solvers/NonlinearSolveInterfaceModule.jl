"""
Interface between data reconciliation problems and [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl):
a Julia package for high-performance and differentiation-enabled nonlinear solvers (Newton methods),
bracketed rootfinding (bisection, Falsi), with sparsity and Newton-Krylov support.
"""
module NonlinearSolveInterfaceModule

using Reexport
using Dictionaries
using LinearAlgebra
using SparseArrays
import NonlinearSolve as NS
using NonlinearSolve: NonlinearFunction, NonlinearLeastSquaresProblem
using NonlinearSolveBase
using ResultTypes: iserror
import SciMLBase

using ..ConstraintsModule
using ..ErrorInfosModule
using ..QCQPModule
using ..ReconciliationModule
using ..PresolveModule: remove_fixed_variables!

@reexport import CommonSolve: init, solve!

export NonlinearSolver

"""
Solve the problem using any solver available from the `NonlinearSolve.jl` Julia library.
"""
struct NonlinearSolver <: ReconciliationSolver
  "Nonlinear solver backend. Use polyalgorithm by default."
  solver::Union{NonlinearSolveBase.AbstractNonlinearSolveAlgorithm, Nothing}
  "Relative tolerance."
  reltol::Float64
  "Absolute tolerance."
  abstol::Float64
  "Maximum number of iterations."
  maxiters::Int64
  "Verbose printing."
  verbose::Bool
end

# Default solver.
function NonlinearSolver(;
  solver = nothing,
  reltol::Float64 = 1e-6,
  abstol::Float64 = 1e-6,
  maxiters::Int64 = 70,
  verbose::Bool = false,
)
  NonlinearSolver(solver, reltol, abstol, maxiters, verbose)
end

# ===================================
# Nonlinear solver pipeline
# ===================================

mutable struct NonlinearPipeline <: ReconciliationPipeline
  "Initial reconciliation problem (not mutated)."
  initial_problem::ReconciliationProblem
  "Reconciliation problem (can be mutated in-place)."
  problem::ReconciliationProblem
  "Solver struct."
  solver::NonlinearSolver
  "Post-process the solution returning infeasible if basic constraints are not satisfied."
  force_basic_constraints::Bool
  "Round results as a post-processing step from the solver."
  round_results::Bool
end

function init(
  problem::ReconciliationProblem,
  solver::NonlinearSolver;
  force_basic_constraints::Bool = true,
  round_results::Bool = false,
)::NonlinearPipeline
  NonlinearPipeline(problem, deepcopy(problem), solver, force_basic_constraints, round_results)
end

"Update the residual term given current values `u`, adding one constraint per equation."
function nlls!(du, u, p, xqp::QuadraticallyConstrainedQuadraticProgram)
  for (i, eq) in enumerate(xqp.constraints)
    du[i] = evaluate_unsafe(eq, u)
  end
end

function jacobian!(J, u, p, xqp::QuadraticallyConstrainedQuadraticProgram)
  for (i, eq) in enumerate(xqp.constraints)
    row = view(J, i, :)
    evaluate_gradient_unsafe!(row, eq, u)
  end
end

function solve!(pipe::NonlinearPipeline)::ReconciliationSolution
  # Formulate the nonlinear problem as a quadratically constrained one.
  (; problem, force_basic_constraints, round_results) = pipe

  # Problem dimensions.
  if !isempty(problem.fixed)
    initial_problem = deepcopy(problem)
    result = remove_fixed_variables!(problem)
    if iserror(result)
      return error_solution(initial_problem, result.error)
    end
  end
  nmeasured = length(problem.measured)
  nunmeasured = length(problem.unmeasured)
  ntot = nmeasured + nunmeasured
  nequations = length(problem.equations)
  @assert ntot > 0 "No variables to solve for."
  @assert nequations > 0 "No equations to solve for."

  # Closure over the quadratic problem.
  xqp = convert(QuadraticallyConstrainedQuadraticProgram, problem)
  f = (du, u, p) -> nlls!(du, u, p, xqp)

  jac = (J, u, p) -> jacobian!(J, u, p, xqp)
  # Unpack solver options.
  (; solver, maxiters, verbose, reltol, abstol) = pipe.solver

  # Needs to be of the same length of the total number of equations. 
  resid_prototype = zeros(Float64, nequations)
  # TODO Evaluate passing the Jacobian explicitly.
  # See: https://docs.sciml.ai/NonlinearSolve/stable/basics/nonlinear_functions/#SciMLBase.NonlinearFunction
  nl_func = NonlinearFunction(f; resid_prototype, jac)
  # TODO Pass custom starting values to the solver.
  u0 = zeros(Float64, ntot)
  nl_prob = NonlinearLeastSquaresProblem(nl_func, u0)

  solver_descr = nothing
  (; value, time) = @timed if isnothing(solver)
    solver_descr = "LevenbergMarquardt + GaussNewton"
    # Stage 1: robust
    lm_maxiters = max(Int64(round(0.65 * maxiters)), 1)
    res = NS.solve(nl_prob, NS.LevenbergMarquardt(); maxiters = lm_maxiters, reltol, abstol)

    # Stage 2: accuracy
    # Call GN if neither of the following holds:
    # 1. Stage 1 succeeded
    # 2. Stage 1 meet strict tolerance
    solver_succeeded = SciMLBase.successful_retcode(res)
    small_residual = maximum(abs, res.resid) < abstol
    if !solver_succeeded && !small_residual
      gn_maxiters = max(maxiters - lm_maxiters, 1)
      NS.solve(nl_prob, NS.GaussNewton(); u0 = res.u, maxiters = gn_maxiters, verbose, reltol, abstol)
    else
      res
    end
  else
    solver_descr = typeof(solver)
    NS.solve(nl_prob, solver; maxiters, verbose, reltol, abstol)
  end
  result = value

  # Fill values from solution. By construction of the QCQP, variables insertion follows the order:
  # measured, unmeasured and fixed (which do not exist at this point since they're removed).
  variables_value = Dictionary{Label, Float64}()
  for (i, k) in enumerate(keys(problem.measured))
    set!(variables_value, k, result.u[i])
  end
  for (i, k) in enumerate(keys(problem.unmeasured))
    set!(variables_value, k, result.u[i + nmeasured])
  end

  status = if result.retcode == NS.ReturnCode.Success
    FEASIBLE
  elseif result.retcode == NS.ReturnCode.MaxIters && (maximum(abs, result.resid) < 1e-1)
    FEASIBLE
  else
    # Status is determined based on analyzing solutions values in the constructor of ReconciliationSolution.
    INFEASIBLE
  end

  # Additional infos.
  infos = Dictionary{String, Any}()
  set!(infos, INFO_NLS_RETCODE, result.retcode)
  set!(infos, INFO_NLS_RESID, result.resid)
  set!(infos, INFO_SV_NLS_TIME, time)
  set!(infos, INFO_SV_SOLVER, solver_descr)

  # Return solution, passing forward solver info.
  ReconciliationSolution(; problem, variables_value, status, infos, force_basic_constraints, round_results)
end

end
