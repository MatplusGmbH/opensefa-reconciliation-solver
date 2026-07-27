"""
Interface between data reconciliation problems and [JuMP.jl](https://github.com/jump-dev/JuMP.jl),
a modeling language for mathematical optimization.

The central type is `JuMPSolver`, which wraps any MOI-compatible optimizer. Solver-specific
convenience constructors are provided in dedicated modules:

- `IpoptSolverModule` — `IpoptSolver()` for smooth nonlinear problems (default).
- `AlpineSolverModule` — `AlpineSolver()` for nonconvex problems requiring a global optimizer.
"""
module JuMPInterfaceModule

using Reexport
using Dictionaries
using LinearAlgebra
using ResultTypes: iserror, unwrap
import JuMP
import MathOptInterface as MOI
import Ipopt

using ..ConstraintsModule
using ..ErrorInfosModule
using ..ReconciliationModule
using ..PresolveModule: find_unobservable_variables_in_constraints, PresolvedReconciliationProblem, remove_fixed_variables!

@reexport import CommonSolve: init, solve!
@reexport import ..ReconciliationModule:
  solver_succeeded, set_starting_values!, is_cost_functional_supported, evaluate_status!

export JuMPSolver,
  add_measured_variables!,
  add_unmeasured_variables!,
  add_fixed_variables!,
  add_equation_constraints!,
  add_objective!,
  build_jump_model!

"""
Holds a solver for a ReconciliationProblem based on the JuMP framework.

The `optimizer` field accepts any MOI-compatible optimizer. The default is `Ipopt.Optimizer`.
Solver-specific convenience constructors (`IpoptSolver`, `AlpineSolver`, ...) are defined in
their respective modules and return a configured `JuMPSolver`.
"""
struct JuMPSolver <: ReconciliationSolver
  "Optimization solver (from MathOptInterface) used."
  optimizer::Union{MOI.AbstractOptimizer, MOI.OptimizerWithAttributes}
  "If 'true', replace fixed variables by the assigned values before solving."
  remove_fixed_variables::Bool
  "If 'true', the problem is treated as feasibilty problem (only solving for constraints) and no data reconciliation is done."
  feasibility::Bool
  "If 'true', constraint equations are relaxed by slack variables that become part of the objective."
  slack_constraints::Bool
  "If 'true', absolute value of unmeasured variable becomes part of the objective."
  slack_unmeasured::Bool
  "Vector of indices for constraint equations that shall be taken into account (no equation IDs but counter within the equation Dictionary). Uses all constraints if 'nothing' (default)."
  subset_constraints::Union{Vector{Int64}, Nothing}
  "Allowed solver tolerance."
  tol::Float64
  "Maximum number of solver iterations."
  max_iter::Int64
  "Dependency detector algorithm eg. mumps."
  dependency_detector::Union{String, Nothing}
  "Use dependency detection with right-hand side (yes, no)."
  dependency_detection_with_rhs::Union{String, Nothing}
  "Compute error propagation and store in the solution infos."
  propagate_errors::Bool
  "Weight factor for quadratic objective (squared norm of unmeasured variables) for feasibility problems."
  objective_regularization_factor::Float64
end

function JuMPSolver(;
  optimizer = Ipopt.Optimizer,
  verbose::Bool = false,
  remove_fixed_variables = true,
  feasibility::Bool = false,
  slack_constraints::Bool = false,
  slack_unmeasured::Bool = false,
  subset_constraints::Union{Vector{Int64}, Nothing} = nothing,
  tol::Float64 = 1e-7,
  max_iter::Int64 = 1000,
  dependency_detector::Union{String, Nothing} = nothing,
  dependency_detection_with_rhs::Union{String, Nothing} = nothing,
  propagate_errors::Bool = true,
  objective_regularization_factor::Float64 = 1e-6,
)
  # If optimizer is already a pre-built MOI.OptimizerWithAttributes (e.g. from AlpineSolver),
  # use it directly to avoid double-wrapping (e.g. adding MOI.Silent() which some optimizers
  # like Alpine do not support).
  moi_optimizer = if optimizer isa MOI.OptimizerWithAttributes
    optimizer
  else
    MOI.OptimizerWithAttributes(optimizer, MOI.Silent() => !verbose)
  end
  JuMPSolver(
    moi_optimizer,
    remove_fixed_variables,
    feasibility,
    slack_constraints,
    slack_unmeasured,
    subset_constraints,
    tol,
    max_iter,
    dependency_detector,
    dependency_detection_with_rhs,
    propagate_errors,
    objective_regularization_factor,
  )
end

function set_starting_values!(initial_values::Dictionary{Label, Float64}, solver::JuMPSolver)
  for (name, val) in pairs(initial_values)
    set!(solver.initial_values, name, val)
  end
end

mutable struct JuMPPipeline{RP <: AbstractReconciliationProblem} <: ReconciliationPipeline
  "Initial problem (not mutated)."
  initial_problem::RP
  "Problem (can be mutated in-place)."
  problem::RP
  "JuMP solver."
  solver::JuMPSolver
  "JuMP model."
  model::Union{Nothing, JuMP.Model}
  "Initial value for each variable."
  initial_values::Dictionary{Label, Float64}
  "Post-process the solution returning infeasible if basic constraints are not satisfied."
  force_basic_constraints::Bool
  "Round results as a post-processing step from the solver."
  round_results::Bool
  "Remove unmeasured variables that are not uniquely defined by the constraints from the solution."
  remove_unobservable_variable_values::Bool
end

"This solver supports additional cost functionals."
is_cost_functional_supported(::JuMPSolver) = true

function init(
  problem::ReconciliationProblem,
  solver::JuMPSolver;
  initial_values::Dictionary{Label, Float64} = Dictionary{Label, Float64}(),
  force_basic_constraints::Bool = true,
  round_results::Bool = true,
  remove_unobservable_variable_values::Bool = false,
)::JuMPPipeline
  JuMPPipeline(
    problem,
    deepcopy(problem),
    solver,
    nothing,
    initial_values,
    force_basic_constraints,
    round_results,
    remove_unobservable_variable_values,
  )
end

"Mapping of variables' names to JuMP variables."
const VariableMap = Dictionary{Label, Pair{Symbol, Int64}}

# Declare the measured variables.
function add_measured_variables!(model::JuMP.Model, rec::ReconciliationProblem, assignment::VariableMap)::JuMP.Model
  nmeasured = length(rec.measured)
  if iszero(nmeasured)
    return model
  end
  JuMP.@variable(model, x[1:nmeasured])
  for (i, n) in enumerate(keys(rec.measured))
    var = get_variable(rec, n)
    # If lower or upper bounds are infinity, it will effectively behave as if no bound was defined.
    lb = lower_bound(var)
    if !ismissing(lb) && !isinf(lb)
      JuMP.set_lower_bound(x[i], lb)
    end
    ub = upper_bound(var)
    if !ismissing(ub) && !isinf(ub)
      JuMP.set_upper_bound(x[i], ub)
    end
    # Populate assignments dictionary.
    set!(assignment, n, :x => i)
  end
  model
end

function add_unmeasured_variables!(
  model::JuMP.Model,
  rec::ReconciliationProblem,
  assignment::VariableMap;
  slack_unmeasured::Bool = false,
)::JuMP.Model
  nunmeasured = length(rec.unmeasured)
  if iszero(nunmeasured)
    return model
  end
  JuMP.@variable(model, y[1:nunmeasured])
  if slack_unmeasured
    JuMP.@variable(model, λu[1:nunmeasured] .>= 1e-3)
  end
  for (i, n) in enumerate(keys(rec.unmeasured))
    var = get_variable(rec, n)
    lb = lower_bound(var)
    if !ismissing(lb) && !isinf(lb)
      JuMP.set_lower_bound(y[i], lb)
    end
    ub = upper_bound(var)
    if !ismissing(ub) && !isinf(ub)
      JuMP.set_upper_bound(y[i], ub)
    end

    if slack_unmeasured
      JuMP.@constraint(model, model[:y][i] <= model[:λu][i])
      JuMP.@constraint(model, model[:y][i] >= -model[:λu][i])
    end

    # Populate assignments dictionary.
    set!(assignment, n, :y => i)
  end
  model
end

"Add fixed variables to the problem (only applicable in case such trivial replacements weren't performed)."
function add_fixed_variables!(model::JuMP.Model, rec::ReconciliationProblem, assignment::VariableMap)::JuMP.Model
  nfixed = length(rec.fixed)
  if iszero(nfixed)
    return model
  end
  JuMP.@variable(model, z[1:nfixed])
  for (i, n) in enumerate(keys(rec.fixed))
    var = get_variable(rec, n)
    lb = lower_bound(var)
    if !ismissing(lb) && !isinf(lb)
      JuMP.set_lower_bound(z[i], lb)
    end
    ub = upper_bound(var)
    if !ismissing(ub) && !isinf(ub)
      JuMP.set_upper_bound(z[i], ub)
    end
    # Populate assignments dictionary.
    set!(assignment, n, :z => i)
  end
  model
end

function add_equation_constraints!(
  model::JuMP.Model,
  rec::ReconciliationProblem,
  assignment::VariableMap;
  slack_constraints::Bool = false,
  subset_constraints::Union{Vector{Int64}, Nothing} = nothing,
  equation_assignment::Dictionary{Int64, Vector{JuMP.ConstraintRef}} = Dictionary{Int64, Vector{JuMP.ConstraintRef}}(),
)::JuMP.Model
  if slack_constraints
    nequations = length(rec.equations)
    JuMP.@variable(model, λc[1:nequations] >= 0.0)
  end

  # For each equation add a constraint, which can be either linear or nonlinear.
  for (i, ei) in enumerate(rec.equations)
    if !isnothing(subset_constraints)
      # Skip all constraints not in the subset index.
      if i ∉ subset_constraints
        @debug "Skipping constraint with index $i"
        continue
      end
    end
    if !iszero(ei.constant_term) && isempty(ei.linear_terms) && isempty(ei.bilinear_terms)
      # Typically by removal of fixed constraints we find an equation with only constants.
      # In any case it shouldn't be added to the model.
      s = evaluate(ei.constant_term)
      if !iszero(s)
        @warn "Equation expected to be zero, got $s."
      end
      continue
    end

    acc = zero(Float64)

    # Accumulate constant terms.
    acc += evaluate(ei.constant_term)

    # Accumulate linear terms.
    for t in ei.linear_terms
      sym, idx = assignment[t.name]
      var = model[sym][idx]
      val = t.factor
      acc += val * var
    end

    # Accumulate bilinear terms.
    for t in ei.bilinear_terms
      sym, idx = assignment[t.name1]
      var1 = model[sym][idx]

      sym, idx = assignment[t.name2]
      var2 = model[sym][idx]

      val = t.factor
      acc += val * var1 * var2
    end

    # Add the constraint.
    # The dictionary `equation_assignment` maps equations of the reconciliation problem to JuMP constraints.
    # Note that the same original equation is related to two equations in the JuMP model, so we store array values.
    ub = slack_constraints ? λc[i] : zero(Float64)
    lb = slack_constraints ? -λc[i] : zero(Float64)

    c = JuMP.@constraint(model, acc <= ub)
    push!(get!(equation_assignment, i, JuMP.ConstraintRef[]), c)
    c = JuMP.@constraint(model, acc >= lb)
    push!(get!(equation_assignment, i, JuMP.ConstraintRef[]), c)
  end

  model
end

function add_initial_values!(
  model::JuMP.Model,
  ::ReconciliationProblem,
  assignment::VariableMap,
  initial_values::Dictionary{Label, Float64},
  variables_value::Dictionary{Label, Float64},
)::JuMP.Model

  # Set a starting values for every variable for which there is an initial value available.
  for (n, val) in pairs(initial_values)
    if haskey(variables_value, n)
      # Variables that have been already assigned can ignore the initial value that was given by the user.
      continue
    end
    if !haskey(assignment, n)
      throw(ArgumentError("The JuMP model is missing variable $n"))
    end
    sym, idx = assignment[n]
    var = model[sym][idx]
    JuMP.set_start_value(var, val)
  end
  model
end

"""
Add the reconciliation optimization objective. If there are no measured variables, set 1.0 as the
objective and the problem should be understood as a feasibility one.

To add linear cost terms on top of the quadratic objective, wrap the problem in an
`OptimizationProblem` (from `OptimizationModule`) and dispatch on that type.
"""
function add_objective!(
  model::JuMP.Model,
  rec::ReconciliationProblem,
  assignment::Dictionary{Label, Pair{Symbol, Int64}};
  feasibility::Bool = false,
  slack_constraints::Bool = false,
  slack_unmeasured::Bool = false,
  objective_regularization_factor = 1e-6,
)::JuMP.Model
  obj = if feasibility || isempty(rec.measured)
    # Constant of one just for the sake of a constant objective.
    y = model[:y]
    objective_regularization_factor * dot(y, y)
  else
    # Weights matrix containing uncertainties.
    Q = weights_matrix(rec)
    # Objective function: minimize (x-xm)^T * Q^{-1} * (x-xm)
    xm = [xi.value for xi in rec.measured]
    x = model[:x]
    dot(x - xm, Q \ (x - xm))
  end
  α = 1.0
  β = 1.0
  if slack_constraints
    obj += α * sum(model[:λc])
  end
  if slack_unmeasured
    obj += β * sum(model[:λu])
  end
  JuMP.@objective(model, Min, obj)
  model
end

"""
Build a JuMP model for a reconciliation problem.

Encapsulates the common model-building logic (variables, constraints, objective) used by
both `ReconciliationProblem` and `OptimizationProblem` workflows.
The objective dispatch is on `full_problem`, so `OptimizationProblem` can override it
by defining a method for `add_objective!(model, ::OptimizationProblem, ...)`.

Returns `(assignment, equation_assignment)`.
"""
function build_jump_model!(
  model::JuMP.Model,
  full_problem::AbstractReconciliationProblem,
  solver::JuMPSolver,
  initial_values::Dictionary{Label, Float64},
  variables_value::Dictionary{Label, Float64},
)::Tuple{VariableMap, Dictionary{Int64, Vector{JuMP.ConstraintRef}}}
  rec = get_reconciliation_problem(full_problem)
  (; slack_unmeasured, slack_constraints, subset_constraints, feasibility, objective_regularization_factor) = solver

  # Add x (measured), y (unmeasured) and z (fixed) variables to the model.
  assignment = VariableMap()
  add_measured_variables!(model, rec, assignment)
  add_unmeasured_variables!(model, rec, assignment; slack_unmeasured)
  add_fixed_variables!(model, rec, assignment)

  # Add initial values (if available).
  if !isempty(initial_values)
    add_initial_values!(model, rec, assignment, initial_values, variables_value)
  end

  # Add equations.
  equation_assignment = Dictionary{Int64, Vector{JuMP.ConstraintRef}}()
  add_equation_constraints!(model, rec, assignment; slack_constraints, subset_constraints, equation_assignment)

  # Add objective function — dispatches on full_problem type.
  add_objective!(
    model,
    full_problem,
    assignment;
    feasibility,
    slack_constraints,
    slack_unmeasured,
    objective_regularization_factor,
  )

  assignment, equation_assignment
end

function extract_values!(
  variables_value::Dictionary{Label, Float64},
  model::JuMP.Model,
  rec::ReconciliationProblem,
)::Dictionary{Label, Float64}
  if !isempty(rec.measured)
    x = model[:x]
    for (i, var) in enumerate(keys(rec.measured))
      set!(variables_value, var, JuMP.value(x[i]))
    end
  end

  if !isempty(rec.unmeasured)
    y = model[:y]
    for (i, var) in enumerate(keys(rec.unmeasured))
      set!(variables_value, var, JuMP.value(y[i]))
    end
  end

  if !isempty(rec.fixed)
    z = model[:z]
    for (i, var) in enumerate(keys(rec.fixed))
      set!(variables_value, var, JuMP.value(z[i]))
    end
  end

  variables_value
end
function extract_values(model::JuMP.Model, rec::ReconciliationProblem)::Dictionary{Label, Float64}
  variables_value = Dictionary{Label, Float64}()
  extract_values!(variables_value, model, rec)
end

function solver_succeeded(sol::ReconciliationSolution, ::JuMPSolver)::Bool
  JuMP.is_solved_and_feasible(sol.infos[INFO_JUMP_MODEL])
end

solver_tolerance(s::JuMPSolver) = s.tol
solver_max_iter(s::JuMPSolver) = s.max_iter

function evaluate_status!(
  jump_status::MOI.TerminationStatusCode;
  infos::Dictionary{String, Any} = Dictionary{String, Any}(),
)::STATUS
  if jump_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED]
    return FEASIBLE
  else
    add_info!(infos, INFO_REASON_INFEASIBILITY, "JuMP Optimization not successful")
    add_info!(infos, INFO_FAILURE_POINT, "JuMP Optimization")
    return INFEASIBLE
  end
end

function solve!(pipe::JuMPPipeline)::ReconciliationSolution
  # Unpack solver settings.
  (;
    problem,
    initial_problem,
    initial_values,
    solver,
    force_basic_constraints,
    round_results,
    remove_unobservable_variable_values,
  ) = pipe
  (; optimizer, tol, max_iter, remove_fixed_variables, dependency_detector, dependency_detection_with_rhs) = solver

  (; value, time) = @timed begin
    # Create optimization instance.
    model = JuMP.Model(optimizer)

    # Update optimizer settings only when supported (e.g. Ipopt supports these; Alpine does not).
    backend = JuMP.backend(model)
    if MOI.supports(backend, MOI.RawOptimizerAttribute("tol"))
      JuMP.set_attribute(model, "tol", tol)
    end
    if MOI.supports(backend, MOI.RawOptimizerAttribute("max_iter"))
      JuMP.set_attribute(model, "max_iter", max_iter)
    end
    if !isnothing(dependency_detector)
      JuMP.set_attribute(model, "dependency_detector", dependency_detector)
    end
    if !isnothing(dependency_detection_with_rhs)
      JuMP.set_attribute(model, "dependency_detection_with_rhs", dependency_detection_with_rhs)
    end

    variables_value = Dictionary{Label, Float64}()
    unobservable_variables = Dictionary{Label, UnmeasuredVariable}()
    unconstrained_variables = Dictionary{Label, MeasuredVariable}()
    rec = get_reconciliation_problem(problem)
    if remove_fixed_variables
      result = remove_fixed_variables!(rec)
      if iserror(result)
        return error_solution(get_reconciliation_problem(initial_problem), result.error)
      end
      prerec = unwrap(result)
      merge!(variables_value, prerec.replacements)
      merge!(unobservable_variables, prerec.unobservable)
      merge!(unconstrained_variables, prerec.unconstrained)
    end

    # Build model.
    assignment, equation_assignment = build_jump_model!(model, problem, solver, initial_values, variables_value)

    # Solve model.
    JuMP.optimize!(model)

    # Extract variable values.
    extract_values!(variables_value, model, rec)

    # Store JuMP's model in the pipeline.
    pipe.model = model

    # Additional infos.
    infos = Dictionary{String, Any}()
    set!(infos, INFO_JUMP_MODEL, model)
    set!(infos, INFO_REDUCED_PROBLEM, rec)
    add_info!(infos, INFO_JUMP_OPTIMIZER, optimizer.optimizer_constructor)
    add_info!(infos, INFO_JUMP_STATUS, JuMP.termination_status(model))
    add_info!(infos, INFO_JUMP_SOLVE_TIME, JuMP.solve_time(model))

    # Set a first status code based on JuMP return.
    status = evaluate_status!(JuMP.termination_status(model); infos)

    # Set variable values for unconstrained measured variables.
    assign_mean!(unconstrained_variables, variables_value)

    # Check for unobservable unmeasured variables.
    merge!(unobservable_variables, find_unobservable_variables_in_constraints(rec, variables_value))
    set!(infos, INFO_N_UNOBSERVABLE, length(unobservable_variables))

    # Optionally remove values of unobservable unmeasured variables.
    if remove_unobservable_variable_values
      delete_unobservable_variable_values!(variables_value, unobservable_variables)
    end
    set!(infos, INFO_UNOBSERVABLE_REMOVED, remove_unobservable_variable_values)

  end
  set!(infos, INFO_SOLVE_TIME, time)

  sol = ReconciliationSolution(;
    problem = get_reconciliation_problem(initial_problem),
    variables_value,
    unobservable_variables,
    status,
    infos,
    force_basic_constraints,
    round_results,
  )

  if JuMP.has_duals(model)
    # Post-process the dual such that its dimension corresponds to the number of equations.
    # The effective dual for a constraint equation 'term = 0' is the sum of the duals for the constraint inequalities 'term <= 0' and 'term >= 0'.
    # The sign changes due to different sign conventions (cf. https://jump.dev/MathOptInterface.jl/stable/background/duality/#Quadratically-Constrained-Quadratic-Programs-(QCQPs) vs. https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions#Nonlinear_optimization_problem.)
    λ = [-(JuMP.dual(c[1]) + JuMP.dual(c[2])) for c in equation_assignment]
    set!(infos, INFO_DUAL, λ)
  end
  if solver.propagate_errors
    uncertainties = propagate_errors_RREF(sol)
    set!(infos, INFO_UNCERTAINTIES, uncertainties)
  end

  sol
end

end
