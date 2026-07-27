"""
Default solution strategy used when calling `solve` without specifying a solver.

The strategy runs three sequential stages:
1. **Presolve** — simplify the problem via structural reductions (see `PresolveModule`).
2. **Starting values** — solve the reduced unmeasured-only system with a nonlinear
   least-squares solver to obtain a good initial point for the optimizer.
3. **Optimization** — solve the full weighted least-squares problem with a nonlinear
   programming solver (default: Ipopt via `JuMPSolver`).

The `DefaultSolver` struct bundles the starting-values solver, the optimization solver,
and a `PresolveConfig`. The `DefaultPipeline` and `DefaultPipelineAfterPresolve` types
track intermediate state across the stages.
"""
module DefaultSolverModule

using Reexport
using Dictionaries
using ResultTypes: @try, iserror, unwrap, Result

using ..ConstraintsModule
using ..ErrorInfosModule
using ..ErrorModule
using ..NonlinearSolveInterfaceModule
using ..JuMPInterfaceModule
using ..PresolveModule
using ..ReconciliationModule
using ..QCQPModule

@reexport import CommonSolve: solve, init, solve!

export DefaultSolver, DefaultPipeline, DefaultPipelineAfterPresolve, PresolveStageResult, StartingReconciliationProblem,
  find_starting_values, check_variables

struct DefaultSolver{S <: Union{ReconciliationSolver, Nothing}, O <: ReconciliationSolver, PC <: PresolveConfig} <:
       ReconciliationSolver
  "Solve to compute starting values to the reconciliation problem."
  starting_values_solver::S
  "Solver for the nonlinear reconciliation optimization problem."
  optimization_solver::O
  "Refers to the presolve config."
  config::PC
end

function DefaultSolver(;
  starting_values_solver::Union{ReconciliationSolver, Nothing} = NonlinearSolver(),
  optimization_solver::ReconciliationSolver = JuMPSolver(; verbose = false, max_iter = 500, propagate_errors = false),
  config::PresolveConfig = PresolveConfig(),
)
  DefaultSolver(starting_values_solver, optimization_solver, config)
end

"If no solver is provided, dispatch on `DefaultSolver`."
function init(problem::ReconciliationProblem; kwargs...)
  init(problem, DefaultSolver(); kwargs...)
end

mutable struct DefaultPipeline <: ReconciliationPipeline
  "Initial reconciliation problem (not mutated)."
  initial_problem::ReconciliationProblem
  "Reconciliation problem (can be mutated in-place)."
  problem::ReconciliationProblem
  "Solver struct."
  solver::DefaultSolver
  "Post-process the solution returning infeasible if basic constraints are not satisfied."
  force_basic_constraints::Bool
  "Round results as a post-processing step from the solver."
  round_results::Bool
  "Remove unmeasured variables that are not uniquely defined by the constraints from the solution."
  remove_unobservable_variable_values::Bool
end

function init(
  problem::ReconciliationProblem,
  solver::DefaultSolver;
  force_basic_constraints::Bool = true,
  remove_unobservable_variable_values::Bool = false,
  round_results::Bool = true,
)::DefaultPipeline
  DefaultPipeline(
    problem,
    deepcopy(problem),
    solver,
    force_basic_constraints,
    round_results,
    remove_unobservable_variable_values,
  )
end

"""
Data structure representing the starting values solution of a reconciliation problem.

Measured variables are replaced by their mean values, therefore the starting values problem
only contains unmeasured variables.
"""
struct StartingReconciliationProblem
  "Presolved system after replacing measured variables. It is nothing if presolve was not applied during this phase."
  prec_unmeasured::Union{PresolvedReconciliationProblem, Nothing}
  "Constraints have been relaxed in order to successfully solve for the unmeasured variables."
  constraints_relaxed::Bool
  "Solution to the presolved system. It is nothing if starting values solver was not applied during this phase, i.e. the problem has been solved before."
  sol_unmeasured::Union{ReconciliationSolution, Nothing}
  "Variable assignments for starting values of the original system."
  starting_values::Dictionary{Label, Float64}
end

"""
Find starting values for the problem, where measured variables are replaced by their mean values.
Unmeasured variables are solved using a nonlinear least squares solver.
"""
function find_starting_values!(
  rec::ReconciliationProblem,
  solver::ReconciliationSolver = NonlinearSolver();
  is_presolved::Bool = false,
  config::PresolveConfig = PresolveConfig(),
  discard_infeasible_equations::Bool = true,
)::Result{StartingReconciliationProblem, SolverError}

  (; time) = @timed begin
    # We assume there are no fixed variables.
    @assert isempty(rec.fixed)

    # Gather all starting values in a single dictionary.
    starting_values = Dictionary{Label, Float64}()

    prec_unmeasured = nothing
    constraints_relaxed = false
    pre = if !isempty(rec.measured) || (isempty(rec.measured) && !is_presolved)
      # Replace each measured variable with its mean value in all equations.
      # Measured variables are replaced by their mean values.
      assign_mean!(rec.measured, starting_values)
      # Idea: measured variables are replaced by its mean values. Then solve the unmeasured variables from the resulting equations.
      # Potential error: mean values don't satisfy the constraint equations
      # Error handling: relax constraints, i.e. remove constraint equations that are not satisfied by the mean values
      # Practical issue: remove_measured_variables! modifies rec even if the removal fails
      # -> relaxing constraints and remove measured variables again can cause conflicts
      # Therefore: Try removal of measured variables with a deepcopy(rec) not modifying the original problem before all conflicting equations are removed.
      rec_for_testing_removal = deepcopy(rec)
      result = remove_measured_variables!(rec_for_testing_removal)
      if iserror(result)
        if discard_infeasible_equations
          while iserror(result) && !isempty(rec.measured)
            # This equation makes the system infeasible, so we remove it.
            eid = result.error.index
            @debug "Discarding infeasible equation $eid."
            remove_equation!(rec, eid)
            constraints_relaxed = true
            rec_for_testing_removal = deepcopy(rec)
            result = remove_measured_variables!(rec_for_testing_removal)
          end
          rec = rec_for_testing_removal
          prec_unmeasured = unwrap(result)
        else
          return result.error
        end
      else
        rec = rec_for_testing_removal
        prec_unmeasured = unwrap(result)
      end

      # Presolve again, eventually simplifying some of the equations.
      @try presolve!(prec_unmeasured; config)

      # Unmeasured variables which were replaced during presolve.
      merge!(starting_values, prec_unmeasured.replacements)
      prec_unmeasured.pre
    else
      # Can skip presolve if the problem has no measured variables and it was presolved before.
      rec
    end

    @assert isempty(pre.fixed)
    @assert isempty(pre.measured)

    # Value nothing indicates no more variables to solve for, and the system is solved at this stage.
    sol_unmeasured = nothing

    if isempty(pre.unmeasured) && !isempty(pre.equations)
      return InternalError("Expected that the starting values system has at least one unmeasured variable.")
    end
  end
  starting_preparation_time = time

  if !isempty(pre.unmeasured)
    (; time) = @timed begin
      # Solve for this new system, which only contains unmeasured variables.
      # Solve using nonlinear least squares root finding for the standing unmeasured variables.
      sol_unmeasured = solve(pre, solver)

      # Starting values solution for remaining unmeasured variables.
      merge!(starting_values, sol_unmeasured.variables_value)
    end

    add_info!(sol_unmeasured.infos, "Starting values preparation time (s)", starting_preparation_time)
    add_info!(sol_unmeasured.infos, "Starting values solving time (s)", time)
  end

  StartingReconciliationProblem(prec_unmeasured, constraints_relaxed, sol_unmeasured, starting_values)
end
function find_starting_values(
  rec::ReconciliationProblem,
  solver::ReconciliationSolver = NonlinearSolver();
  is_presolved::Bool = false,
  config::PresolveConfig = PresolveConfig(),
  discard_infeasible_equations::Bool = true,
)::Result{StartingReconciliationProblem, SolverError}
  find_starting_values!(deepcopy(rec), solver; is_presolved, config, discard_infeasible_equations)
end

# ================
# Presolve stage

"Collects the outputs of the presolve stage."
struct PresolveStageResult
  presolve_data::PresolvedReconciliationProblem
  "Accumulated variable assignments: presolve replacements plus mean values of unconstrained measured variables."
  variables_value::Dictionary{Label, Float64}
  unobservable_variables::Dictionary{Label, UnmeasuredVariable}
end

"If no solver is provided, dispatch on `DefaultSolver`."
function init(problem::ReconciliationProblem, prerec::PresolvedReconciliationProblem; kwargs...)
  variables_value = Dictionary{Label, Float64}()
  merge!(variables_value, prerec.replacements)
  unobservable_variables = prerec.unobservable

  # Set variable values for unconstrained measured variables.
  assign_mean!(prerec.unconstrained, variables_value)

  ps = PresolveStageResult(prerec, variables_value, unobservable_variables)
  init(problem, ps, DefaultSolver(); kwargs...)
end

mutable struct DefaultPipelineAfterPresolve <: ReconciliationPipeline
  "Initial reconciliation problem (not mutated)."
  initial_problem::ReconciliationProblem
  "Reconciliation problem (can be mutated in-place)."
  presolve_result::PresolveStageResult
  "Solver struct."
  solver::DefaultSolver
  "Post-process the solution returning infeasible if basic constraints are not satisfied."
  force_basic_constraints::Bool
  "Round results as a post-processing step from the solver."
  round_results::Bool
  "Remove unmeasured variables that are not uniquely defined by the constraints from the solution."
  remove_unobservable_variable_values::Bool
  "Additional information from the solver backend, including timing results."
  infos::Dictionary{String, Any}
end

function init(
  problem::ReconciliationProblem,
  presolve_result::PresolveStageResult,
  solver::DefaultSolver;
  force_basic_constraints::Bool = true,
  remove_unobservable_variable_values::Bool = false,
  round_results::Bool = true,
  infos::Dictionary{String, Any} = Dictionary{String, Any}(),
)::DefaultPipelineAfterPresolve
  DefaultPipelineAfterPresolve(
    problem,
    presolve_result,
    solver,
    force_basic_constraints,
    round_results,
    remove_unobservable_variable_values,
    infos,
  )
end

"""
Run the presolve stage.

Returns a `PresolveStageResult` on success, or a `ReconciliationSolution` if the problem
is already fully solved at this stage or the presolver fails.
"""
function run_presolve_stage!(
  pipe::DefaultPipeline,
  infos::Dictionary{String, Any},
)::Union{PresolveStageResult, ReconciliationSolution}
  (; initial_problem, problem, solver, force_basic_constraints, round_results) = pipe
  (; config) = solver

  # It is safe to apply the in-place version since the user problem has been copied (`initial_problem`).
  (; value, time) = @timed presolve!(problem; config)
  set!(infos, INFO_PRESOLVE_TIME, time) # TODO Define timings struct.
  if iserror(value)
    set!(infos, INFO_HEURISTIC, "Presolve")
    return error_solution(initial_problem, value.error; infos)
  end
  presolve_data = unwrap(value)

  variables_value = Dictionary{Label, Float64}()
  merge!(variables_value, presolve_data.replacements)
  unobservable_variables = presolve_data.unobservable

  # Set variable values for unconstrained measured variables.
  assign_mean!(presolve_data.unconstrained, variables_value)

  if force_basic_constraints && !has_basic_constraints(problem, variables_value, unobservable_variables)
    # If the basic constraints are not satisfied, always return an infeasibile status.
    add_info!(infos, INFO_REASON_INFEASIBILITY, "Basic constraints violated")
    add_info!(infos, INFO_FAILURE_POINT, "Postprocessing presolve solution")
    return ReconciliationSolution(;
      problem = initial_problem,
      variables_value,
      unobservable_variables,
      status = CONTRADICTORY,
      infos,
      force_basic_constraints,
      round_results,
    )
  end

  # Return if there are no more variables to solve for.
  # This will happen in trivial systems that can be solved already at the presolve stage.
  if isempty(presolve_data.variables_kept)
    @debug "Problem solved at the presolve stage."
    # Sanity checks: all variables from the initial problem need to be either
    # assigned due to presolve steps or otherwise found to be unobservable.
    @assert all(keys(initial_problem.variables)) do var
      (var in presolve_data.variables_removed) || (var in keys(presolve_data.unconstrained)) ||
        (var in keys(presolve_data.unobservable))
    end
    # Check that there are no more variables left in the presolved problem.
    @assert isempty(presolve_data.pre.variables)
    set!(infos, INFO_HEURISTIC, "Presolve")
    return ReconciliationSolution(;
      problem = initial_problem,
      variables_value,
      unobservable_variables,
      status = FEASIBLE,
      infos,
      force_basic_constraints,
      round_results,
      is_global = true,
    )
  end

  PresolveStageResult(presolve_data, variables_value, unobservable_variables)
end

# ======================
# Starting values stage

"""
Collects the outputs of the starting values stage.

The starting values problem is the presolved problem (`PresolveStageResult.presolve_data.pre`)
with measured variables substituted by their mean values. Only unmeasured variables remain,
and they are solved using a nonlinear least squares solver.
"""
struct StartingValuesStageResult
  "Whether the starting values solution should be used directly as the final solution."
  use_solution::Bool
  "Starting values available to seed the optimizer, regardless of whether the starting values solver succeeded."
  starting_values::Union{Dictionary{Label, Float64}, Nothing}
  "Solution for unmeasured variables from the starting values solver, or `nothing` if unavailable."
  sol_unmeasured::Union{ReconciliationSolution, Nothing}
  "Whether the starting values solver converged successfully."
  starting_values_succeeded::Bool
  "Unobservable unmeasured variables identified during the removal of measured variables in `find_starting_values!`."
  unobservable_variables::Dictionary{Label, UnmeasuredVariable}
end

"""
Run the starting values stage.

Computes starting values for unmeasured variables and decides whether they can be used
directly as the final solution.
"""
function run_starting_values_stage(
  pipe::DefaultPipelineAfterPresolve,
  ps::PresolveStageResult,
  infos::Dictionary{String, Any},
)::StartingValuesStageResult
  (; solver, force_basic_constraints) = pipe
  (; config) = solver
  (; presolve_data) = ps

  alg_start = solver.starting_values_solver
  starting = if !isnothing(alg_start)
    (; value, time) = @timed find_starting_values(presolve_data.pre, alg_start; is_presolved = true, config)
    set!(infos, INFO_SV_TIME, time)
    value
  else
    nothing
  end

  # If the starting values solver failed, proceed to the optimizer without making use of starting values.
  starting_values_succeeded = if iserror(starting) || isnothing(starting)
    false
  else
    starting = unwrap(starting)
    # If `sol_unmeasured` does not exist, this indicates that finding the starting values was done only with presolving.
    isnothing(starting.sol_unmeasured) || starting.sol_unmeasured.status == FEASIBLE
  end
  if !isnothing(alg_start)
    set!(infos, INFO_SV_SUCCEEDED, starting_values_succeeded)
  else
    set!(infos, INFO_SV_SUCCEEDED, "Starting values skipped")
  end
  if !iserror(starting) && !isnothing(starting) && !isnothing(starting.sol_unmeasured)
    # Forward timing and solver-identity keys from the starting values sub-solution
    # into the pipeline infos so that the final ReconciliationSolution exposes them
    # at the top level, without callers needing to dig into sol_unmeasured.
    merge_info!(infos, starting.sol_unmeasured.infos, [
      INFO_SV_SOLVER,
      INFO_SV_PREP_TIME,
      INFO_SV_SOLVE_TIME,
      INFO_SV_NLS_TIME,
    ])
  end
  if iserror(starting)
    for key in keys(error_infos(starting.error))
      set!(infos, "Starting values solver - " * key, error_infos(starting.error)[key])
    end
  end

  starting_values = nothing
  sol_unmeasured = nothing
  constraints_relaxed = false
  sv_unobservable_variables = Dictionary{Label, UnmeasuredVariable}()
  if isa(starting, StartingReconciliationProblem)
    # Always capture starting values: used directly when succeeded, or as a seed for the optimizer retry heuristic otherwise.
    starting_values = starting.starting_values
    if starting_values_succeeded
      (; sol_unmeasured, constraints_relaxed) = starting
    end
    # Carry forward unobservable variables identified when measured variables were removed.
    if !isnothing(starting.prec_unmeasured)
      merge!(sv_unobservable_variables, starting.prec_unmeasured.unobservable)
    end
  end

  use_solution = if !starting_values_succeeded || constraints_relaxed
    false
  else
    # If the starting values succeeded and the problem was not relaxed, check whether it is a promising solution.
    starting_status = evaluate_status!(evaluate_score(presolve_data.pre, starting_values))
    # Starting values solution is used only if the solution is good (i.e. status in SUCCESS_CODES)
    # and every variable is solved (i.e. has assigned value or detected as unobservable).
    # Hence PARTIALLY_FEASIBLE is not sufficient to stop computing the solution.
    starting_status in [FEASIBLE, FEASIBLE_WITH_UNOBSERVABLE_VARIABLES]
  end

  # If the starting values shall be used, check if they satisfy the basic constraints (if required).
  if use_solution && force_basic_constraints &&
     !has_basic_constraints(presolve_data.pre, starting_values, presolve_data.unobservable)
    @debug "Starting values don't satisfy the basic constraints."
    use_solution = false
  end

  StartingValuesStageResult(
    use_solution,
    starting_values,
    sol_unmeasured,
    starting_values_succeeded,
    sv_unobservable_variables,
  )
end

# ================
# Optimizer stage

"""
Run the optimizer stage.

Calls the nonlinear programming solver, with an optional retry using starting values.
Returns the updated `variables_value`, `unobservable_variables`, and `status`.

`variables_value_base` must contain the presolve replacements so that `evaluate_score`
can be computed against the full initial problem.
"""
function run_optimizer_stage(
  pipe::DefaultPipelineAfterPresolve,
  ps::PresolveStageResult,
  variables_value_base::Dictionary{Label, Float64},
  infos::Dictionary{String, Any};
  sv::StartingValuesStageResult = StartingValuesStageResult(
    false,
    nothing,
    nothing,
    false,
    Dictionary{Label, UnmeasuredVariable}(),
  ),
)
  (; initial_problem, solver, remove_unobservable_variable_values, round_results) = pipe
  alg_nl = solver.optimization_solver

  # TODO Add algorithm used to the solution struct and remove from infos.
  tol = solver_tolerance(alg_nl)
  isnothing(tol) || set!(infos, INFO_SOLVER_TOLERANCE, tol)
  max_iter = solver_max_iter(alg_nl)
  isnothing(max_iter) || set!(infos, INFO_MAX_ITERATIONS, max_iter)

  (; value, time) = if sv.starting_values_succeeded
    @timed solve(
      ps.presolve_data.pre,
      alg_nl;
      initial_values = sv.starting_values,
      remove_unobservable_variable_values,
      round_results,
    )
  else
    @timed solve(ps.presolve_data.pre, alg_nl; remove_unobservable_variable_values, round_results)
  end
  set!(infos, INFO_SOLVE_TIME, time)
  sol = value
  if haskey(sol.infos, INFO_JUMP_SOLVE_TIME)
    add_info!(infos, INFO_JUMP_SOLVE_TIME, sol.infos[INFO_JUMP_SOLVE_TIME])
  end
  optimization_solver_succeeded = solver_succeeded(sol, alg_nl)

  # Heuristic: if optimization solver did not succeed but starting values are available
  # and haven't been used yet, try using them. See `Muehl_2024_2` for an example.
  if !sv.starting_values_succeeded && !optimization_solver_succeeded && !isnothing(sv.starting_values)
    @debug "Solving second time"
    (; value, time) =
      @timed solve(
        ps.presolve_data.pre,
        alg_nl;
        initial_values = sv.starting_values,
        remove_unobservable_variable_values,
        round_results,
      )
    set!(infos, INFO_SOLVE2_TIME, time)
    sol = value
    if haskey(sol.infos, INFO_JUMP_SOLVE_TIME)
      add_info!(infos, INFO_JUMP_SOLVE2_TIME, sol.infos[INFO_JUMP_SOLVE_TIME])
    end
    optimization_solver_succeeded = solver_succeeded(sol, alg_nl)
  end

  # Accumulate diagnostic keys from the optimizer solution.
  # Use merge_info! (add_info! semantics) so that values from the retry attempt,
  # if it occurs, stack rather than silently overwrite earlier diagnostics.
  merge_info!(
    infos,
    sol.infos,
    [
      INFO_JUMP_OPTIMIZER,
      INFO_JUMP_STATUS,
      INFO_N_UNOBSERVABLE,
      INFO_UNOBSERVABLE_REMOVED,
      INFO_ERROR_MESSAGE,
      INFO_REASON_INFEASIBILITY,
      INFO_PROBLEMATIC_REPLACEMENTS,
      INFO_CONFLICTING_REPLACEMENTS,
      INFO_CONFLICTING_EQ_DEPS,
      INFO_FAILURE_POINT,
      INFO_PREVIOUS_REPLACEMENTS,
    ],
  )
  # Overwrite (not accumulate) structural results — these are single-valued by design.
  # "Dual" and "Uncertainties" may have been partially written by an earlier stage;
  # the optimizer result is authoritative and must replace, not append.
  overwrite_info!(infos, sol.infos, [INFO_REDUCED_PROBLEM, INFO_UNCERTAINTIES, INFO_DUAL])

  unobservable_variables = sol.unobservable_variables
  variables_value = copy(variables_value_base)
  # The computed values are passed even if the low-level `optimization_solver_succeeded` returns `false`
  # and the status is evaluated later
  merge!(variables_value, sol.variables_value)

  # Heuristic: if the final QCQP solve was successful then use it; otherwise
  # stick to the starting solve one (with measured variables at their mean values).
  if optimization_solver_succeeded
    @debug "Using optimization solver values."
    set!(infos, INFO_HEURISTIC, "Optimizer")
  end

  (; status = sol.status, variables_value, unobservable_variables)
end

# ===============
# Default solver

function solve!(pipe::DefaultPipeline)::ReconciliationSolution
  (; initial_problem, solver, force_basic_constraints, round_results, remove_unobservable_variable_values) = pipe
  infos = Dictionary{String, Any}()

  # Stage 1: presolve.
  ps = run_presolve_stage!(pipe, infos)
  ps isa ReconciliationSolution && return ps

  # Continue with Stage 2
  pipe_after_presolve = DefaultPipelineAfterPresolve(
    initial_problem,
    ps,
    solver,
    force_basic_constraints,
    round_results,
    remove_unobservable_variable_values,
    infos,
  )
  solve!(pipe_after_presolve)
end

function solve!(pipe::DefaultPipelineAfterPresolve)::ReconciliationSolution
  (; initial_problem, presolve_result, force_basic_constraints, round_results, remove_unobservable_variable_values, infos) =
    pipe
  ps = presolve_result

  # Return if there are no more variables to solve for.
  if isempty(ps.presolve_data.pre.variables)
    set!(infos, INFO_HEURISTIC, "Presolve")
    return ReconciliationSolution(;
      problem = initial_problem,
      variables_value = ps.variables_value,
      unobservable_variables = ps.unobservable_variables,
      status = FEASIBLE,
      infos,
      force_basic_constraints,
      round_results,
    )
  end

  # Stage 2: starting values.
  sv = run_starting_values_stage(pipe, ps, infos)

  variables_value = ps.variables_value
  unobservable_variables = ps.unobservable_variables

  if sv.use_solution
    @debug "Problem solved in the starting values stage."
    set!(infos, INFO_HEURISTIC, "Starting values")
    merge!(variables_value, sv.starting_values)
    # Check for unobservable unmeasured variables (usually done after solving with JuMP).
    # Also include unobservable variables identified during removal of measured variables in find_starting_values!.
    merge!(unobservable_variables, sv.unobservable_variables)
    merge!(unobservable_variables, find_unobservable_variables_in_constraints(ps.presolve_data, variables_value))
    set!(infos, INFO_N_UNOBSERVABLE, length(unobservable_variables))
    # Optionally remove unobservable unmeasured variables.
    if remove_unobservable_variable_values
      delete_unobservable_variable_values!(variables_value, unobservable_variables)
    end
    set!(infos, INFO_UNOBSERVABLE_REMOVED, remove_unobservable_variable_values)
    status = if isnothing(sv.sol_unmeasured)
      FEASIBLE
    else
      @assert sv.sol_unmeasured.status in [FEASIBLE, FEASIBLE_WITH_UNOBSERVABLE_VARIABLES] # otherwise the solution shouldn't be used
      if isempty(unobservable_variables)
        FEASIBLE
      else
        FEASIBLE_WITH_UNOBSERVABLE_VARIABLES
      end
    end
  else
    # Stage 3: optimizer.
    opt = run_optimizer_stage(pipe, ps, variables_value, infos; sv)
    status = if !isempty(unobservable_variables) && (opt.status == FEASIBLE)
      FEASIBLE_WITH_UNOBSERVABLE_VARIABLES
    else
      opt.status
    end
    variables_value = opt.variables_value
    merge!(unobservable_variables, opt.unobservable_variables)
  end

  # Add presolved problem to infos if not done before (e.g. for starting values solution).
  if !haskey(infos, INFO_REDUCED_PROBLEM)
    set!(infos, INFO_REDUCED_PROBLEM, ps.presolve_data.pre)
  end

  finalsol =
    ReconciliationSolution(;
      problem = initial_problem,
      variables_value,
      unobservable_variables,
      status,
      infos,
      force_basic_constraints,
      round_results,
    )

  # Error propagation if not done before (e.g. for starting values solution).
  if finalsol.status in SUCCESS_CODES && !haskey(finalsol.infos, INFO_UNCERTAINTIES)
    uncertainties = propagate_errors_RREF(finalsol)
    set!(finalsol.infos, INFO_UNCERTAINTIES, uncertainties)
  end

  finalsol
end

end
