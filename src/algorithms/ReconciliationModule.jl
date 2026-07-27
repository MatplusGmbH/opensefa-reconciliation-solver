"""
Core types and algorithms for data reconciliation.

Data reconciliation finds the most likely values of measured quantities subject to model
constraints by minimizing a weighted sum of squared adjustments (a least-squares problem).

Key types:
- `ReconciliationProblem` — holds variables (measured, unmeasured, fixed) and constraint equations.
- `ReconciliationSolution` — result of a solve call, including variable values and solution status.
- `ReconciliationSolver` — abstract type for solver backends.
- `ReconciliationPipeline` — abstract type for multi-stage solution pipelines.
"""
module ReconciliationModule

using Reexport
using Dictionaries
using LinearAlgebra
using ResultTypes: @try, iserror, unwrap, Result
using SparseArrays
using StructEquality

using ..ConstraintsModule
using ..STANModule
using ..ErrorModule
using ..ComparisonModule
using ..RowEchelonModule
using ..SparseRowEchelonModule
using ..ErrorInfosModule

@reexport import CommonSolve: init, solve!
@reexport import ..ConstraintsModule: is_linear, is_bilinear, evaluate, is_transfer_coefficient
@reexport import ..ErrorInfosModule: add_info!

export
  AbstractReconciliationProblem,
  ReconciliationProblem,
  ReconciliationSolution,
  ReconciliationSolver,
  ReconciliationPipeline,
  get_reconciliation_problem,
  SolutionScore,
  SolutionComparison,
  get_variable,
  add_variable!,
  add_equation!,
  remove_variable!,
  remove_equation!,
  add_linear_term!,
  add_bilinear_term!,
  update_vars_dicts!,
  error_solution,
  weights_matrix,
  evaluate_objective,
  evaluate_fixed_discrepancy,
  evaluate_score,
  evaluate_status!,
  evaluate_globality,
  has_basic_constraints,
  solver_succeeded,
  solver_tolerance,
  solver_max_iter,
  set_starting_values!,
  is_cost_functional_supported,
  propagate_errors,
  propagate_errors_RREF,
  jacobians_and_hessian,
  jacobians,
  compare,
  jacobian,
  STATUS,
  SUCCESS_CODES,
  PARTIALLY_FEASIBLE,
  FEASIBLE,
  FEASIBLE_WITH_UNOBSERVABLE_VARIABLES,
  INFEASIBLE,
  CONTRADICTORY,
  UNKNOWN,
  DEFAULT_ZTOL,
  EID

const DEFAULT_ZTOL = default_tolerance(Float64).ztol

# ===================
# Problem definition
# ===================

"Index for a given equation. Indices are unique within the same system."
const EID = Int64

"""
Abstract type for reconciliation problems.

Allows different problem types (e.g., ReconciliationProblem, OptimizationProblem)
to share common solver infrastructure while maintaining type safety.
"""
abstract type AbstractReconciliationProblem end

"""
Contains various dictionaries that hold problem data for efficient lookups.

To create an empty problem, use `ReconciliationProblem()`.
The following interface functions should be used to act upon reconciliation problems
(instead of operations on the arrays directly):

- `add_variable!` -- Add a variable to the problem.
- `remove_variable!` -- Remove a variable from the problem. 
- `add_equation!` -- Add an equation to the problem.
- `remove_equation!` -- Remove an equation from the problem.

Each of these methods works for an individual variable or an array.
"""
@kwdef struct ReconciliationProblem <: AbstractReconciliationProblem
  "Model's name."
  name::Label = ""
  "Vector of equations."
  equations::Dictionary{EID, Equation} = Dictionary{EID, Equation}()
  "Mapping variable names to their corresponding type."
  variables::Dictionary{Label, Type{<:Variable}} = Dictionary{Label, Type{<:Variable}}()
  "Vector of measured variables (known with uncertainty)."
  measured::Dictionary{Label, MeasuredVariable} = Dictionary{Label, MeasuredVariable}()
  "Vector of unmeasured variables (unknowns or model parameters)."
  unmeasured::Dictionary{Label, UnmeasuredVariable} = Dictionary{Label, UnmeasuredVariable}()
  "Vector of variables which for which fixed values have been assigned (known without uncertainty)."
  fixed::Dictionary{Label, FixedVariable} = Dictionary{Label, FixedVariable}()
  "Mapping of variable names to the indices of equations in which they participate."
  var_to_eq::Dictionary{Label, Vector{EID}} = Dictionary{Label, Vector{EID}}() # initialize_variable_map(variables, equations)
  "Reverse lookup for each equation to the variables that participate in it."
  eq_to_var::Dictionary{EID, Vector{Label}} = Dictionary{EID, Vector{Label}}()
end

function ReconciliationProblem(model::String)
  convert(ReconciliationProblem, parse(STANTrace, model))
end

function Base.show(io::IO, rec::ReconciliationProblem)
  name = rec.name
  nequations = length(rec.equations)
  nunmeasured = length(rec.unmeasured)
  nmeasured = length(rec.measured)
  nfixed = length(rec.fixed)
  nvariables = nunmeasured + nmeasured + nfixed

  println(io, "• Data reconciliation optimization problem with:")
  println(io, "  • Name: $name")
  println(io, "  • Equations: $nequations")
  println(io, "  • Variables: $nvariables")
  println(io, "    • Measured: $nmeasured")
  println(io, "    • Unmeasured: $nunmeasured")
  println(io, "    • Fixed: $nfixed")
end

@struct_equal ReconciliationProblem

"Each variable is mapped to the indices of the equations in which it participates (ie. it is in a linear, bilinear term or both)."
function initialize_variable_map(
  variables::Dictionary{Label, Type{<:Variable}},
  equations::Dictionary{EID, Equation},
)::Dictionary{Label, Vector{EID}}
  var_to_eq = Dictionary{Label, Vector{EID}}()
  label_buffer = Label[]
  for (i, eq) in pairs(equations)
    empty!(label_buffer)
    for var in term_names!(label_buffer, eq)
      push!(get!(var_to_eq, var, Int64[]), i)
    end
  end
  var_to_eq
end

"""
    get_reconciliation_problem(problem::AbstractReconciliationProblem)

Extract the underlying `ReconciliationProblem` from any `AbstractReconciliationProblem`.
For `ReconciliationProblem` itself, returns self. For wrapper types like `OptimizationProblem`,
this must be overloaded to extract the wrapped `ReconciliationProblem`.
"""
get_reconciliation_problem(rec::ReconciliationProblem) = rec

# =========================
# Problem utility methods
# =========================

"Return if the system only has linear flow equations."
function is_linear(rec::ReconciliationProblem)::Bool
  all(is_linear, rec.equations)
end

"Return if the system has at least one bilinear equation."
function is_bilinear(rec::ReconciliationProblem)::Bool
  !is_linear(rec)
end

"Return if the named variable represents a transfer coefficient."
function is_transfer_coefficient(rec::ReconciliationProblem, name::Label)::Bool
  if !haskey(rec.variables, name)
    throw(ArgumentError("Variable $name does not exist in the problem."))
  end
  var = get_variable(rec, name)
  is_transfer_coefficient(var)
end

"""
Add an equation to the reconciliation problem.

If `strict` is `true` (default) an error is thrown unless all of the following conditions hold:

- All variables in the equation need to exist in the problem (see `add_variable!`), otherwise an error is thrown.
- All linear terms have distinct names.
- All bilinear terms have distinct names, modulo squared terms and permutations within the same term.

Such assumptions are useful to improve the efficiency of certain algorithms.
"""
function add_equation!(rec::ReconciliationProblem, eq::Equation; strict::Bool = true)::ReconciliationProblem
  names = term_names(eq)
  if strict
    # Ensure that all variables in the equation are already defined in the reconciliation problem.
    all_variables_exist = names ⊆ keys(rec.variables)
    if !all_variables_exist
      throw(ArgumentError("The equation contains variables that are not defined in the reconciliation problem."))
    end

    # Ensure that there are no different linear terms with the same name.
    found = Set{Label}()
    all_lt_distinct = true
    for lt in eq.linear_terms
      if lt.name in found
        throw(ArgumentError("The equation contains repeated linear terms, reduce them first."))
        break
      else
        push!(found, lt.name)
      end
    end
    if !all_lt_distinct
      throw(ArgumentError("The equation contains repeated linear terms, reduce them first."))
    end

    # Ensure that there are no different bilinear terms with the same name, modulo squared terms and permutations.
    found = Set{Tuple{Label, Label}}()
    all_bt_distinct = true
    for bt in eq.bilinear_terms
      if ((bt.name1, bt.name2) in found) || ((bt.name2, bt.name1) in found)
        all_bt_distinct = false
        break
      else
        push!(found, (bt.name1, bt.name2))
      end
    end
    if !all_bt_distinct
      throw(ArgumentError("The equation contains repeated bilinear terms, reduce them first."))
    end
  end

  # Add the equation to the problem.
  m = length(rec.equations) + 1
  insert!(rec.equations, m, eq)

  # For this equation, update variables mapping that appear on each term.
  for var in names
    push!(get!(rec.var_to_eq, var, EID[]), m)
  end

  # Update reverse lookup from equations to variables.
  insert!(rec.eq_to_var, m, names)

  rec
end
function add_equation!(rec::ReconciliationProblem, eqs::Vector{Equation}; strict::Bool = true)::ReconciliationProblem
  foreach(eqs) do eq
    add_equation!(rec, eq; strict)
  end
  rec
end

"""
Remove an equation from the reconciliation problem.
The method does not perform any variable removals.
"""
function remove_equation!(rec::ReconciliationProblem, id::EID)::ReconciliationProblem
  delete!(rec.equations, id)
  for var in rec.eq_to_var[id]
    # For each variable appearing in the equation, update the variable to EID map. 
    pos = findfirst(==(id), rec.var_to_eq[var])
    @assert !isnothing(pos)
    deleteat!(rec.var_to_eq[var], pos)
  end
  delete!(rec.eq_to_var, id)
  rec
end
function remove_equation!(rec::ReconciliationProblem, indices::Vector{EID})::ReconciliationProblem
  foreach(indices) do id
    remove_equation!(rec, id)
  end
  rec
end

"Add a new variable to the reconciliation problem."
function add_variable!(rec::ReconciliationProblem, var::MeasuredVariable)::ReconciliationProblem
  # Using `insert!` instead of `set!` implies that an error is thrown if the variable already exists.
  insert!(rec.variables, var.name, MeasuredVariable)
  insert!(rec.measured, var.name, var)
  insert!(rec.var_to_eq, var.name, EID[])
  rec
end
function add_variable!(rec::ReconciliationProblem, var::UnmeasuredVariable)::ReconciliationProblem
  insert!(rec.variables, var.name, UnmeasuredVariable)
  insert!(rec.unmeasured, var.name, var)
  insert!(rec.var_to_eq, var.name, EID[])
  rec
end
function add_variable!(rec::ReconciliationProblem, var::FixedVariable)::ReconciliationProblem
  insert!(rec.variables, var.name, FixedVariable)
  insert!(rec.fixed, var.name, var)
  insert!(rec.var_to_eq, var.name, EID[])
  rec
end
function add_variable!(rec::ReconciliationProblem, vars::Vector{<:Variable})::ReconciliationProblem
  foreach(vars) do var
    add_variable!(rec, var)
  end
  rec
end

"""
Remove a variable from the reconciliation problem.

Note: removing a variable does not affect the equations in which it may participate; those need to be
adjusted accordingly.

In some cases we use this method to workaround the fact that variable structs are immutable. For example,
in order to change the bound of a variable we can add it and remove one with the same name. By doing so,
since the name is the same the coupling with the equations in which it participates should remain unchanged.
"""
function remove_variable!(rec::ReconciliationProblem, name::Label)::ReconciliationProblem
  # The method does not modify `eq_to_var`, ie. it is assumed that the variable does not appear in any equation.
  T = get(rec.variables, name, nothing)
  if isnothing(T)
    throw(ArgumentError("Variable $name doesn't exist."))
  end
  # Remove from general variable dictionaries, then from the specific one.
  delete!(rec.variables, name)
  if haskey(rec.var_to_eq, name)
    delete!(rec.var_to_eq, name)
  end
  if T == MeasuredVariable
    delete!(rec.measured, name)
  elseif T == UnmeasuredVariable
    delete!(rec.unmeasured, name)
  elseif T == FixedVariable
    delete!(rec.fixed, name)
  else
    throw(ArgumentError("Variable type undefined."))
  end
  rec
end
function remove_variable!(rec::ReconciliationProblem, names::Vector{Label})::ReconciliationProblem
  foreach(names) do name
    remove_variable!(rec, name)
  end
  rec
end

"Retrieve the variable instance from a name, or nothing if the variable is not in the problem."
function get_variable(rec::ReconciliationProblem, name::Label)::Union{Variable, Nothing}
  T = get(rec.variables, name, nothing)
  isnothing(T) && return
  if T == MeasuredVariable
    get(rec.measured, name, nothing)
  elseif T == UnmeasuredVariable
    get(rec.unmeasured, name, nothing)
  elseif T == FixedVariable
    get(rec.fixed, name, nothing)
  else
    throw(ArgumentError("Variable type undefined."))
  end
end

"""
Add a linear term 'factor * var' to the equation 'eid' of a ReconciliationProblem 'rec'
assuming that there is no linear term with variable 'var' in equation 'eid' yet.
"""
function add_new_linear_term!(rec::ReconciliationProblem, eid::EID, var::Label, factor::Float64)
  eq = rec.equations[eid]
  ltnew = LinearTerm(var, factor)
  push!(eq.linear_terms, ltnew)
  # Update var_to_eq and eq_to_var if the variable didn't participate in the equation before
  if !(var in rec.eq_to_var[eid])
    push!(rec.eq_to_var[eid], var)
    push!(rec.var_to_eq[var], eid)
  end
end

"""
Add a bilinear term 'factor * var1 * var2' to the equation 'eid' of a ReconciliationProblem 'rec'
assuming that there is no bilinear term with variables 'var1' and 'var2' in equation 'eid' yet.
"""
function add_new_bilinear_term!(rec::ReconciliationProblem, eid::EID, var1::Label, var2::Label, factor::Float64)
  eq = rec.equations[eid]
  btnew = BilinearTerm(var1, var2, factor)
  push!(eq.bilinear_terms, btnew)
  # Update var_to_eq and eq_to_var if the variables didn't participate in the equation before
  if !(var1 in rec.eq_to_var[eid])
    push!(rec.eq_to_var[eid], var1)
    push!(rec.var_to_eq[var1], eid)
  end
  if !(var2 in rec.eq_to_var[eid])
    push!(rec.eq_to_var[eid], var2)
    push!(rec.var_to_eq[var2], eid)
  end
end

"""
Add 'value' to the factor of the LinearTerm with index 'lt_index' of the equation 'eid' of a ReconciliationProblem 'rec'.
If the added 'value' and the previous factor cancel each other, the LinearTerm is removed from the equation and the variable stored in 'removed_variables'.
Note that the dictionaries rec.var_to_eq and rec.eq_to_var and rec.variables are not updated. (Requires further checking using 'removed_variables'.)
"""
function add_to_linear_term!(
  rec::ReconciliationProblem,
  eid::EID,
  lt_index::Int64,
  value::Float64,
  removed_variables::Vector{Label} = Label[],
)
  eq = rec.equations[eid]
  previous_linear_term = eq.linear_terms[lt_index]
  var = previous_linear_term.name
  previous_factor = previous_linear_term.factor
  if !isapproxzero(value + previous_factor)
    eq.linear_terms[lt_index] = LinearTerm(var, value + previous_factor)
  else
    # The existing and new linear terms cancel each other
    deleteat!(eq.linear_terms, lt_index)
    if var ∉ removed_variables
      push!(removed_variables, var)
    end
  end
end

"""
Add 'value' to the factor of the BilinearTerm with index 'bt_index' of the equation 'eid' of a ReconciliationProblem 'rec'.
If the added 'value' and the previous factor cancel each other, the BilinearTerm is removed from the equation and the variables stored in 'removed_variables'.
Note that the dictionaries rec.var_to_eq and rec.eq_to_var and rec.variables are not updated. (Requires further checking using 'removed_variables'.)
"""
function add_to_bilinear_term!(
  rec::ReconciliationProblem,
  eid::EID,
  bt_index::Int64,
  value::Float64,
  removed_variables::Vector{Label} = Label[],
)
  eq = rec.equations[eid]
  previous_bilinear_term = eq.bilinear_terms[bt_index]
  var1 = previous_bilinear_term.name1
  var2 = previous_bilinear_term.name2
  previous_factor = previous_bilinear_term.factor
  if !isapproxzero(value + previous_factor)
    eq.bilinear_terms[bt_index] = BilinearTerm(var1, var2, value + previous_factor)
  else
    # The existing and new linear terms cancel each other
    deleteat!(eq.bilinear_terms, bt_index)
    push!(removed_variables, var1, var2)
  end
end

"""
Add a linear term 'factor * var' to the equation 'eid' of a ReconciliationProblem 'rec'.
If a linear term with variable 'var' existed before but cancels with the added linear term, the variable 'var' is stored in 'removed_variables'.
"""
function add_linear_term!(
  rec::ReconciliationProblem,
  eid::EID,
  var::Label,
  factor::Float64,
  removed_variables::Vector{Label} = Label[],
)
  @assert haskey(rec.equations, eid)
  @assert haskey(rec.variables, var)
  eq = rec.equations[eid]
  # Check if a term with varible 'var' already exists in a linear term in this equations
  idx_existing = get_index(eq.linear_terms, var)
  if isnothing(idx_existing)
    # A linear term with this name does not exist, so we add a new one.
    add_new_linear_term!(rec, eid, var, factor)
  else
    # A linear term with this name already exists, so we add the new value to the previous one.
    add_to_linear_term!(rec, eid, idx_existing, factor, removed_variables)
  end
end

"""
Add a biinear term 'factor * var1 * var2' to the equation 'eid' of a ReconciliationProblem 'rec'.
If a bilinear term with variables 'var1' and 'var2' existed before but cancels with the added bilinear term, the variables 'var1' and 'var2' are stored in 'removed_variables'.
"""
function add_bilinear_term!(
  rec::ReconciliationProblem,
  eid::EID,
  var1::Label,
  var2::Label,
  factor::Float64,
  removed_variables::Vector{Label} = Label[],
)
  @assert haskey(rec.equations, eid)
  @assert haskey(rec.variables, var1)
  @assert haskey(rec.variables, var2)
  eq = rec.equations[eid]
  # Check if 'var1 * var2' or 'var2 * var1' already exists in a bilinear term in this equations
  idx_existing = get_index(eq.bilinear_terms, var1, var2)
  if isnothing(idx_existing)
    # A bilinear term with these variables does not exist, so we add a new one.
    add_new_bilinear_term!(rec, eid, var1, var2, factor)
  else
    # A bilinear term with these variables already exists, so we add the new value to the previous one.
    add_to_bilinear_term!(rec, eid, idx_existing, factor, removed_variables)
  end
end

"""
Check for variables in 'vars_to_check' and equations in 'eqs_to_check' if they have to be removed from var_to_eq and eq_to_var in the ReconciliationProblem 'rec'.
If necessary, remove from the dictionaries and store removed variables in 'removed_variables'.
Note that it is not checked if a variable is unconstrained/free. (Requires further checking using 'removed_variables'.)
If 'vars_to_check' or 'eqs_to_check' is not specified (i.e., 'nothing'), check all variables and equations from the ReconciliationProblem.
"""
function update_vars_dicts!(
  rec::ReconciliationProblem,
  vars_to_check::Union{Vector{Label}, Nothing} = nothing,
  eqs_to_check::Union{Vector{EID}, Nothing} = nothing,
  removed_variables::Vector{Label} = Label[],
)
  if isnothing(vars_to_check)
    vars_to_check = collect(keys(rec.variables))
  end
  if isnothing(eqs_to_check)
    eqs_to_check = collect(keys(rec.equations))
  end
  # Remove variables from var_to_eq and eq_to_var if they don't participate anymore
  for varname in vars_to_check
    for eid in eqs_to_check
      eq = rec.equations[eid]
      if !(varname in term_names(eq)) && varname in rec.eq_to_var[eid]
        deleteat!(rec.eq_to_var[eid], findfirst(==(varname), rec.eq_to_var[eid]))
        deleteat!(rec.var_to_eq[varname], findfirst(==(eid), rec.var_to_eq[varname]))
        push!(removed_variables, varname)
      end
    end
  end
end

"Convert a `STANTrace` into a reconciliation problem."
function Base.convert(::Type{ReconciliationProblem}, trace::STANTrace)
  rec = ReconciliationProblem(; name = trace.name)
  for var in trace.initial_variables
    add_variable!(rec, var)
  end
  for eq in trace.equations
    add_equation!(rec, eq; strict = true)
  end
  rec
end

"""
Evaluate the fixed constraints (fixed variables of the reconciliation problem).
Returns the element-wise absolute difference between expected (fixed) values and results. 
"""
function evaluate_fixed_discrepancy(rec::ReconciliationProblem, vars::Vector{<:Variable})::Dictionary{Label, Float64}
  result = Dictionary{Label, Float64}()
  for (name, var) in pairs(rec.fixed)
    # Map each variable from result into its corresponding fixed variable.
    idx = findfirst(x -> x.name == var.name, vars)
    if isnothing(idx)
      throw(ArgumentError("Measured variable $name not found in the result."))
    end
    # Store the element-wise absolute difference.
    set!(result, var.name, abs(var.value - vars[idx].value))
  end
  result
end
function evaluate_fixed_discrepancy(
  rec::ReconciliationProblem,
  vars::Dictionary{Label, Float64},
)::Dictionary{Label, Float64}
  result = Dictionary{Label, Float64}()
  for (name, var) in pairs(rec.fixed)
    # Store the element-wise absolute difference between the expected value and the computed value.
    set!(result, var.name, abs(var.value - vars[name]))
  end
  result
end

"""
Evaluate the equations of the reconciliation problem.
Returns the numerical evaluation of the left-hand side of each equation (it is expected to be zero).
"""
function evaluate(rec::ReconciliationProblem, results::Dictionary{Label, Float64})::Vector{Float64}
  lhs = fill(zero(Float64), length(rec.equations))
  for (i, eq) in enumerate(rec.equations)
    lhs[i] = evaluate(eq, results)
  end
  lhs
end

"""
Return the weights matrix `Q` containing variance-covariance components of the measured quantities.
If the problem has no measured variables, `nothing` is returned.
"""
function weights_matrix(rec::ReconciliationProblem)::Union{Diagonal, Nothing}
  if isempty(rec.measured)
    return
  end
  uncertainties = [x.uncertainty^2 for x in rec.measured]
  Diagonal(uncertainties)
end

"""
Evalutes the weighted least-squares function `(x - xm)^T * Qinv * (x - xm)` between measured and reconciled variables,
where `Q` is the matrix of variance-covariance uncertainties.
If the problem has no measured variables, `nothing` is returned.
"""
function evaluate_objective(rec::ReconciliationProblem, vars::Vector{<:Variable})::Union{Float64, Nothing}
  # Compute the weights matrix and the quadratic function.
  Q = weights_matrix(rec)
  if isnothing(Q)
    return
  end
  Qinv = inv(Q)

  # Compute `y = x - xm`.
  y = Vector{Float64}()
  for (name, var) in pairs(rec.measured)
    idx = findfirst(x -> x.name == var.name, vars)
    if isnothing(idx)
      throw(ArgumentError("Measured variable $(x.name) not found in the result."))
    end
    push!(y, vars[idx].value - var.value)
  end
  dot(y, Qinv, y)
end
function evaluate_objective(rec::ReconciliationProblem, vars::Dictionary{Label, Float64})::Union{Float64, Nothing}
  # Compute the weights matrix and the quadratic function.
  Q = weights_matrix(rec)
  if isnothing(Q)
    return
  end

  # Compute `y = x - xm`.
  y = Vector{Float64}()
  for (name, var) in pairs(rec.measured)
    if !haskey(vars, name)
      throw(ArgumentError("Measured variable $name not found in the result."))
    end
    push!(y, vars[name] - var.value)
  end
  dot(y, Q \ y)
end

# =================
# Solver pipeline
# =================

"Extend this type for solver backends."
abstract type ReconciliationSolver end

"Return the solver tolerance for logging, or `nothing` if not applicable."
solver_tolerance(::ReconciliationSolver) = nothing

"Return the maximum iteration count for logging, or `nothing` if not applicable."
solver_max_iter(::ReconciliationSolver) = nothing

"Extend this type for solver pipelines."
abstract type ReconciliationPipeline end

"""
Status returned by the solver.

- `FEASIBLE`                           : All variables have assigned value (no unobservable variables)
- `PARTIALLY_FEASIBLE`                 : Some variables have assigned value (from calling Presolve)
- `FEASIBLE_WITH_UNOBSERVABLE_VARIABLES` : There are some variables which are unobservable
- `CONTRADICTORY`                      : Cases for which we find infeasibility after variable replacements
- `INFEASIBLE`                         : Cases in which the solver was not able to find the solution
- `UNKNOWN`                            : Default value for methods (never returned by a method)
"""
@enum STATUS UNKNOWN FEASIBLE PARTIALLY_FEASIBLE FEASIBLE_WITH_UNOBSERVABLE_VARIABLES CONTRADICTORY INFEASIBLE

Base.convert(::Type{String}, x::STATUS) = string(x)

# Default implementation.
function set_starting_values!(initial_values::Dictionary{Label, Float64}, sol::ReconciliationSolver)
  throw(ArgumentError("Not implemented."))
end

"Any of the following is considered a successful solve."
const SUCCESS_CODES = (FEASIBLE, PARTIALLY_FEASIBLE, FEASIBLE_WITH_UNOBSERVABLE_VARIABLES)


"""
Holds scores that describe the quality of a ReconciliationSolution.
"""
struct SolutionScore
  "Percentage of variables which got a value assigned or identified as unobservable."
  solved_variables_score::Float64
  "Percentage of variables for which a value was assigned."
  assigned_variables_score::Float64
  "Percentage of variables for which a unique solution value was found (ignoring unobservable variables)."
  solved_variables_unique_score::Float64
  "Percentage of equations which can be evaluated with the solution values."
  solved_equations_score::Float64
  "Dictionary holding the evaluation for the equations after variable replacements. (Equations aren't evaluated if variables missing a value occur therein.)"
  equations_value::Dictionary{EID, Float64}
  "Dictionary holding the absolute residual for each equation."
  absolute_residuals::Dictionary{EID, Float64}
  "The residual (infinity norm of absolute_residuals) or nothing if the problem has no equations."
  residual::Union{Float64, Nothing}
  "Dictionary holding the relative residual for each equation."
  relative_residuals::Dictionary{EID, Float64}
  "The maximum relative residual (infinity norm of relative_residuals) or nothing if the problem has no equations."
  max_relative_residual::Union{Float64, Nothing}
end

function Base.show(io::IO, score::SolutionScore)
  solved_variables = round(score.solved_variables_score, sigdigits = 5)
  assigned_variables = round(score.assigned_variables_score, sigdigits = 5)
  solved_variables_unique = round(score.solved_variables_unique_score, sigdigits = 5)
  solved_equations = round(score.solved_equations_score, sigdigits = 5)
  residual = if !isnothing(score.residual)
    round(score.residual, sigdigits = 5)
  else
    nothing
  end
  residual_str = isnothing(residual) ? "N/A" : string(residual)

  max_relative_residual = if !isnothing(score.max_relative_residual)
    round(score.max_relative_residual, sigdigits = 5)
  else
    nothing
  end
  max_relative_residual_str = isnothing(max_relative_residual) ? "N/A" : string(max_relative_residual)

  println(io, "• Solution score:")
  println(io, "    • Percentage of solved variables: $solved_variables")
  println(io, "    • Percentage of assigned variables: $assigned_variables")
  println(io, "    • Percentage of uniquely solved variables: $solved_variables_unique")
  println(io, "    • Percentage of solved equations: $solved_equations")
  println(io, "    • Residual: $residual_str")
  println(io, "    • Max relative residual: $max_relative_residual_str")
end


"""
Compute the relative residual for a single equation.

The relative residual is `|r(x)| / (|c| + ||q|| ||x|| + ||P|| ||x||^2)` where:
- `r(x)` is the residual (signed equation value)
- `c` is the constant term
- `q` are the linear term coefficients
- `P` are the bilinear term coefficients

This computation relies on `equations_value`, i.e., it can only be computed if the equation
can be evaluated. This typically corresponds to all variables having assigned values.
As an exception to such rule, for a bilinear term 'factor * var1 * var2' we allow that a variable is unassigned if the other is zero.

Note: The constraint equation does not define P uniquely but only P_{ij} + P_{ji} for indices i, j.
We always set P such that for i != j it holds either P_{ij} == 0 or P_{ji} == 0, but theoretically
it could be done in a different way (e.g., defining P to be symmetric).

The formula uses ||x|| computed from ALL unique variables appearing in the equation.
This ensures consistency: the same variable vector is used for both linear and bilinear contributions.

This function computes the denominator directly from the equation structure without QCQP conversion.
"""
function compute_relative_residual(
  eq::Equation,
  residual::Float64,
  variables_value::Dictionary{Label, Float64},
)::Float64
  # Constant term contribution: |c|.
  c_norm = abs(eq.constant_term.value)

  # Collect ALL unique variable names from both linear and bilinear terms.
  all_var_names = term_names(eq)

  # Compute ||x|| for ALL variables in the equation.
  x_norm_sq = zero(Float64)
  for var_name in all_var_names
    if haskey(variables_value, var_name)
      x_norm_sq += variables_value[var_name]^2
    end
  end
  x_norm = sqrt(x_norm_sq)

  # Linear term contribution: ||q|| ||x||.
  q_norm_sq = zero(Float64)
  for t in eq.linear_terms
    q_norm_sq += t.factor^2
  end
  q_norm = sqrt(q_norm_sq)
  linear_contribution = q_norm * x_norm

  # Bilinear term contribution: ||P|| ||x||^2.
  P_norm_sq = zero(Float64)
  for t in eq.bilinear_terms
    @assert !iszero(t.factor) "Bilinear term factor cannot be zero as the equation has been evaluated"
    P_norm_sq += t.factor^2
  end
  P_norm = sqrt(P_norm_sq)
  bilinear_contribution = P_norm * x_norm^2

  # Compute denominator.
  denominator = c_norm + linear_contribution + bilinear_contribution

  # Avoid division by zero.
  if iszero(denominator)
    # If denominator is zero, return absolute residual.
    return abs(residual)
  end

  abs(residual) / denominator
end

"""
Given the initial problem and a dictionary of assigned values (eg. result of optimization), then
this function estimates the percentage of variables (and equations) whose value is assigned.
"""
function evaluate_score(
  rec::ReconciliationProblem,
  variables_value::Dictionary{Label, Float64},
  unobservable::Dictionary{Label, UnmeasuredVariable} = Dictionary{Label, UnmeasuredVariable}(),
)::SolutionScore
  nvars_calculated =
    length([var for var in keys(rec.variables) if var in keys(variables_value) || var in keys(unobservable)])
  nvars_assigned = length(variables_value)
  nvars_assigned_unique = length(setdiff(keys(variables_value), keys(unobservable)))
  nvars = length(rec.variables)
  @assert !iszero(nvars)
  solved_variables_score = nvars_calculated / nvars
  assigned_variables_score = nvars_assigned / nvars
  solved_variables_unique_score = nvars_assigned_unique / nvars

  # Evaluate each equation assiging the value after variable replacements.
  equations_value = Dictionary{EID, Float64}()

  for (i, eq) in pairs(rec.equations)
    # If the equation's evaluation fails due to missing variable values, such equation is skept affecting `solved_equations_score`.
    try
      val = evaluate(eq, variables_value)
      set!(equations_value, i, val)
    catch
      @debug "Equation $i has missing variables from the results vector, skipping."
    end
  end
  neqs = length(rec.equations)
  neqs_calculated = length(equations_value)
  solved_equations_score = iszero(neqs) ? 1.0 : neqs_calculated / neqs

  # Compute the absolute residual for each equation.
  absolute_residuals = map(abs, equations_value)

  residual = if isempty(absolute_residuals)
    nothing
  else
    maximum(absolute_residuals)
  end

  # Compute the relative residual for each equation.
  relative_residuals = Dictionary{EID, Float64}()
  for (i, val) in pairs(equations_value)
    eq = rec.equations[i]
    rel_res = compute_relative_residual(eq, val, variables_value)
    set!(relative_residuals, i, rel_res)
  end

  max_relative_residual = if isempty(relative_residuals)
    nothing
  else
    maximum(abs, relative_residuals)
  end

  SolutionScore(
    solved_variables_score,
    assigned_variables_score,
    solved_variables_unique_score,
    solved_equations_score,
    equations_value,
    absolute_residuals,
    residual,
    relative_residuals,
    max_relative_residual,
  )
end

"Checks for a given SolutionScore if every equation satisfies either absolute_residual < atol or relative_residual < rtol."
function check_residuals(score::SolutionScore; atol::Float64, rtol::Float64)::Bool
  check = true
  for (eid, absres) in pairs(score.absolute_residuals)
    if absres >= atol && abs(score.relative_residuals[eid]) >= rtol
      check = false
      break
    end
  end
  check
end

function evaluate_status!(
  score::SolutionScore;
  status::STATUS = UNKNOWN,
  infos::Dictionary{String, Any} = Dictionary{String, Any}(),
)::STATUS
  if status in [UNKNOWN, FEASIBLE, FEASIBLE_WITH_UNOBSERVABLE_VARIABLES, PARTIALLY_FEASIBLE]
    atol = default_tolerance(Float64).ztol
    rtol = default_tolerance(Float64).rtol
    if check_residuals(score; atol, rtol) == false
      status = INFEASIBLE
      add_info!(infos, INFO_REASON_INFEASIBILITY, "Residual too big")
      add_info!(infos, INFO_FAILURE_POINT, "Processing solver return")
    elseif score.solved_variables_unique_score == 1
      status = FEASIBLE
    elseif score.solved_variables_score == 1
      status = FEASIBLE_WITH_UNOBSERVABLE_VARIABLES
    else
      status = PARTIALLY_FEASIBLE
    end
  end
  return status
end

"Evaluates if it can be proven that any local solution of a given ReconciliationProblem is global."
function evaluate_globality(rec::ReconciliationProblem, status::STATUS)::Bool
  if status in SUCCESS_CODES
    if isempty(rec.measured)
      # Problem without measured variables (feasibility problem)
      return true
    elseif is_linear(rec)
      # All constraint equations are linear, hence the problem is convex
      return true
    end
  end
  # Always return false if no proper solution was found
  return false
end

"""
Holds the reconciliation problem as well as values after variable reconciliation.
Such values can be used to judge the quality of the solution.
"""
struct ReconciliationSolution
  "Reconciliation problem (includes also problems without measured variables, so feasibility and not reconciliation strictly speaking)."
  problem::ReconciliationProblem
  "Value of the objective or nothing if the system has no measured variables. It is missing if it cannot be evaluated due to missing variable assignments."
  objective_value::Union{Float64, Nothing, Missing}
  "Assignment after of values to each variable."
  variables_value::Dictionary{Label, Float64}
  "Unobservable variables (i.e. variables that cannot be assigned uniquely)"
  unobservable_variables::Dictionary{Label, UnmeasuredVariable}
  "Solver status."
  status::STATUS
  "Information on the quality of the solution."
  score::SolutionScore
  "Indicates if it can be proven that the solution is global."
  is_global::Bool
  "Additional information from the solver backend, including timing results."
  infos::Dictionary{String, Any}
end

function ReconciliationSolution(;
  problem::ReconciliationProblem,
  variables_value::Dictionary{Label, Float64},
  unobservable_variables::Dictionary{Label, UnmeasuredVariable} = Dictionary{Label, UnmeasuredVariable}(),
  status::STATUS = UNKNOWN,
  infos::Dictionary{String, Any} = Dictionary{String, Any}(),
  force_basic_constraints::Bool = true,
  round_results::Bool = false,
  is_global::Union{Bool, Nothing} = nothing,
)
  if round_results
    # Redefine any variable that saturates its lower (or upper) bound.
    for name in keys(variables_value)
      var = get_variable(problem, name)
      val = variables_value[name]
      lb, ub = lower_bound(var), upper_bound(var)
      if !ismissing(lb) && relaxed_isapprox(lb, val) && !tight_leq(lb, val)
        set!(variables_value, name, lb)
      end
      if !ismissing(ub) && relaxed_isapprox(ub, val) && !tight_geq(ub, val)
        set!(variables_value, name, ub)
      end
    end
  end

  # Evaluate solution score.
  sol_score = evaluate_score(problem, variables_value, unobservable_variables)

  # Additional information.
  can_evaluate_objective = all(v -> v in keys(variables_value), keys(problem.measured))
  objective_value = if can_evaluate_objective
    evaluate_objective(problem, variables_value)
  else
    missing
  end
  set!(infos, INFO_COMPUTED_VARIABLES, length(variables_value))
  set!(infos, INFO_MISSING_VARIABLES, length(problem.variables) - length(variables_value))
  set!(infos, INFO_COMPUTED_EQUATIONS, length(sol_score.equations_value))
  set!(infos, INFO_MISSING_EQUATIONS, length(problem.equations) - length(sol_score.equations_value))
  set!(infos, INFO_COMPUTED_VARIABLES_EXTREMA, extrema(variables_value, init = zeros(2)))
  tc_value = [v for (k, v) in pairs(variables_value) if is_transfer_coefficient(problem, k)]
  set!(infos, INFO_TC_EXTREMA, isempty(tc_value) ? nothing : extrema(tc_value))

  # Use reduced problem for subsequent checkings if available
  rec_reduced = haskey(infos, INFO_REDUCED_PROBLEM) ? infos[INFO_REDUCED_PROBLEM] : problem
  # Don't rely on the JuMP status for feasibility problems
  if haskey(infos, INFO_JUMP_STATUS)
    status = evaluate_status!(infos[INFO_JUMP_STATUS]; infos)
  end
  if isempty(rec_reduced.measured) && !isequal(status, CONTRADICTORY)
    status = UNKNOWN
  end

  if force_basic_constraints
    # If the basic constraints are not satisfied, always return an infeasibile status.
    if !has_basic_constraints(problem, variables_value, unobservable_variables)
      if !isequal(status, CONTRADICTORY)
        status = INFEASIBLE
      end
      add_info!(infos, INFO_REASON_INFEASIBILITY, "Basic constraints violated")
      add_info!(infos, INFO_FAILURE_POINT, "Postprocessing solution")
    end
  end

  # Update solution status if the solver didn't explictly pass it.
  status = evaluate_status!(sol_score; status, infos)

  # Update globality information if the solver didn't explicitly pass it.
  if isnothing(is_global)
    is_global = evaluate_globality(rec_reduced, status)
  end


  ReconciliationSolution(
    problem,
    objective_value,
    variables_value,
    unobservable_variables,
    status,
    sol_score,
    is_global,
    infos,
  )
end

function Base.show(io::IO, rsol::ReconciliationSolution)
  name = rsol.problem.name
  nequations = length(values(rsol.problem.equations))
  nunmeasured = length(values(rsol.problem.unmeasured))
  nmeasured = length(values(rsol.problem.measured))
  nfixed = length(values(rsol.problem.fixed))
  nvariables = nunmeasured + nmeasured + nfixed
  obj = rsol.objective_value
  if !isnothing(obj)
    obj = round(rsol.objective_value, sigdigits = 4)
  end
  (; status, is_global, score) = rsol

  println(io, "• Solution of data reconciliation optimization problem with:")
  println(io, "  • Name: $name")
  println(io, "  • Equations: $nequations")
  println(io, "  • Variables: $nvariables")
  println(io, "    • Measured: $nmeasured")
  println(io, "    • Unmeasured: $nunmeasured")
  println(io, "    • Fixed: $nfixed")
  println(io, "  • Objective value: $obj")
  println(io, "  • Status: $status")
  println(io, "  • Globality proven: $is_global")
  println(io, "  $score")
end

"Add an info to a ReconciliationSolution. Check if the info key already exists."
function add_info!(rsol::ReconciliationSolution, key::String, new_info::Any)
  add_info!(rsol.infos, key, new_info)
end

function Base.convert(
  ::Type{ReconciliationSolution},
  trace::STANTrace;
  use_initial_results::Bool = false,
  use_kelly::Bool = false,
)
  problem = convert(ReconciliationProblem, trace)

  # Evaluate the objective function.
  results = if use_initial_results
    trace.initial_results
  elseif use_kelly
    trace.kelly_results
  else
    trace.final_results
  end

  # Evaluate each variable according to the results array.
  variables_value = results.variables_value
  unobservable_variables = results.unassigned_variables
  infos = Dictionary{String, Any}()
  set!(infos, INFO_HEURISTIC, "STAN solution")
  merge_info!(infos, trace.infos)
  set!(infos, INFO_N_UNOBSERVABLE, length(unobservable_variables))

  # Determine the solution status.
  if !isempty(results.unassigned_variables)
    set!(infos, INFO_UNOBSERVABLE_REMOVED, true)
    @debug "There exist variables without assigned value in the results vector."
  else
    set!(infos, INFO_UNOBSERVABLE_REMOVED, false)
  end

  set!(infos, INFO_UNCERTAINTIES, results.uncertainties)

  ReconciliationSolution(; problem, variables_value, unobservable_variables, infos)
end

function error_solution(
  problem::ReconciliationProblem,
  err::Union{InfeasibleEquationError, InternalError};
  infos::Dictionary{String, Any} = Dictionary{String, Any}(),
)::ReconciliationSolution
  variables_value = Dictionary{Label, Float64}()
  if err isa ContradictionError
    status = CONTRADICTORY
  else
    status = INFEASIBLE
  end
  add_info!(infos, INFO_ERROR_MESSAGE, err.msg)
  merge_info!(infos, err.infos)
  ReconciliationSolution(; problem, variables_value, status, infos)
end

# ==================
# Property checking
# ==================

function evaluate_fixed_discrepancy(rsol::ReconciliationSolution)::Dictionary{Label, Float64}
  vars = [FixedVariable(var, val) for (var, val) in pairs(rsol.variables_value)]
  evaluate_fixed_discrepancy(rsol.problem, vars)
end

"Extend this method for specific solvers."
function solve!(pipe::ReconciliationPipeline)::ReconciliationSolution
  throw(ArgumentError("Not implemented for $(typeof(pipe))"))
end

"Whether the solver succeeded to produce an answer (total or partial)."
function solver_succeeded(sol::ReconciliationSolution, ::ReconciliationSolver)::Bool
  # Default implementation.
  sol.status in SUCCESS_CODES
end

"""
Check if the solution satisfies the basic constraints:

- All measured and unmeasured variables have assigned nonnegative values.
- All transfer coefficients have assigned a value in the `[0, 1]` interval.
- Upper and lower bounds for each variable are respected.
"""
function has_basic_constraints(sol::ReconciliationSolution)::Bool
  problem = sol.problem
  has_basic_constraints(problem, sol.variables_value, sol.unobservable_variables)
end
function has_basic_constraints(
  rec::ReconciliationProblem,
  variables_value::Dictionary{Label, Float64},
  unobservable::Dictionary{Label, UnmeasuredVariable} = Dictionary{Label, UnmeasuredVariable}(),
)::Bool
  result = true
  for var_name in keys(rec.variables)
    var = get_variable(rec, var_name)

    if !haskey(variables_value, var_name)
      if haskey(unobservable, var_name)
        # If the variable is unobservable then do not check that it is present in the variables assignment.
        @debug "The variable $(var_name) from the original problem is unobservable, no value assigned."
      else
        @debug "The variable $(var_name) from the original problem does not have an assigned value."
      end
      continue
    end

    val = variables_value[var_name]
    # Check that lower and upper bounds are respected.
    result = satisfies_bounds(var, val)
    if result == false
      break
    end
  end
  result
end

"By default assume that cost functionals are not supported."
function is_cost_functional_supported(::ReconciliationSolver)
  false
end


"""
Holds the result of the comparision of two ReconciliationSolutions.
"""
struct SolutionComparison
  "Solutions that have been compared."
  solutions::Tuple{ReconciliationSolution, ReconciliationSolution}
  "Name of the underlying ReconciliationProblem."
  problem_name::String
  "Variables that have assigned values in both ReconciliationSolutions."
  shared_variables::Dictionary{Label, Type{<:Variable}}
  "Absolute deviation of variable values."
  vars_absolute_deviation::Dictionary{Label, Float64}
  "Relative deviation of variable values."
  vars_relative_deviation::Dictionary{Label, Float64}
  "Indicates if variables are identified as observable/not observable by both solutions."
  vars_equal_observability::Dictionary{Label, Bool}
  "Equations that have assigned values in both ReconciliationSolutions."
  shared_equations::Dictionary{EID, Equation}
  "Absolute deviation of equation values."
  eqs_absolute_deviation::Dictionary{EID, Float64}
  "Tolerance in which absolute deviation of variables is neglected."
  abstol::Float64
end

function SolutionComparison(
  solutions::Tuple{ReconciliationSolution, ReconciliationSolution},
  problem_name::String,
  shared_variables::Dictionary{Label, Type{<:Variable}},
  vars_absolute_deviation::Dictionary{Label, Float64},
  vars_relative_deviation::Dictionary{Label, Float64},
  vars_equal_observability::Dictionary{Label, Bool},
  shared_equations::Dictionary{EID, Equation},
  eqs_absolute_deviation::Dictionary{EID, Float64};
  abstol::Float64 = 1.0e-7,
)
  SolutionComparison(
    solutions,
    problem_name,
    shared_variables,
    vars_absolute_deviation,
    vars_relative_deviation,
    vars_equal_observability,
    shared_equations,
    eqs_absolute_deviation,
    abstol,
  )
end

function Base.show(io::IO, solcomp::SolutionComparison)
  probname = solcomp.problem_name
  (sol1, sol2) = solcomp.solutions
  nvars1 = length(sol1.problem.variables)
  nvars2 = length(sol2.problem.variables)
  nvars_shared = length(solcomp.shared_variables)
  nvars_equal_observability = count(solcomp.vars_equal_observability)
  neqs1 = length(sol1.problem.equations)
  neqs2 = length(sol2.problem.equations)
  neqs_shared = length(solcomp.shared_equations)
  if !iszero(nvars_shared)
    max_var_dev_abs = round(maximum(solcomp.vars_absolute_deviation), sigdigits = 5)
    max_var_dev_rel = round(maximum(solcomp.vars_relative_deviation), sigdigits = 5)
  else
    max_var_dev_abs = "-"
    max_var_dev_rel = "-"
  end
  vars_is_unique = dictionary(
    collect(
      var => (!haskey(sol1.unobservable_variables, var) && !haskey(sol2.unobservable_variables, var)) for
      var in keys(solcomp.shared_variables)
    ),
  )
  vars_observable_absolute_deviation = dictionary(
    collect(
      var => solcomp.vars_absolute_deviation[var] for var in keys(solcomp.shared_variables) if vars_is_unique[var]
    ),
  )
  vars_observable_relative_deviation = dictionary(
    collect(
      var => solcomp.vars_relative_deviation[var] for var in keys(solcomp.shared_variables) if vars_is_unique[var]
    ),
  )
  if !iszero(count(vars_is_unique))
    max_observable_var_dev_abs = round(maximum(vars_observable_absolute_deviation), sigdigits = 5)
    max_observable_var_dev_rel = round(maximum(vars_observable_relative_deviation), sigdigits = 5)
  else
    max_observable_var_dev_abs = "-"
    max_observable_var_dev_rel = "-"
  end
  vars_relative_deviation_essential = dictionary(
    collect(
      var => solcomp.vars_relative_deviation[var] for
      var in keys(solcomp.vars_absolute_deviation) if solcomp.vars_absolute_deviation[var] > solcomp.abstol
    ),
  )
  vars_observable_relative_deviation_essential = dictionary(
    collect(
      var => solcomp.vars_relative_deviation[var] for
      var in keys(vars_relative_deviation_essential) if vars_is_unique[var]
    ),
  )
  if !iszero(length(vars_relative_deviation_essential))
    max_var_dev_rel_ess = round(maximum(vars_relative_deviation_essential), sigdigits = 5)
  else
    max_var_dev_rel_ess = 0.0
  end
  if !iszero(length(vars_observable_relative_deviation_essential))
    max_observable_var_dev_rel_ess = round(maximum(vars_observable_relative_deviation_essential), sigdigits = 5)
  else
    max_observable_var_dev_rel_ess = 0.0
  end
  if !iszero(neqs_shared)
    max_eq_dev_abs = round(maximum(solcomp.eqs_absolute_deviation), sigdigits = 5)
  else
    max_eq_dev_abs = "-"
  end

  println(io, "• Comparision of two ReconciliationSolutions:")
  println(io, "  • Underlying ReconciliationProblem: $probname")
  println(io, "  • Variables...")
  if isequal(nvars1, nvars2)
    println(io, "    • ... in underlying problem: $nvars1")
  else
    println(io, "    • ... in underlying problem for the first solution: $nvars1")
    println(io, "    • ... in underlying problem for the second solution: $nvars2")
  end
  println(io, "    • ... shared by both solutions: $nvars_shared")
  println(io, "    • ... with equal observability identification: $nvars_equal_observability")
  println(io, "    • Maximum absolute deviation: $max_var_dev_abs")
  println(io, "    • Maximum absolute deviation among observable variables: $max_observable_var_dev_abs")
  println(io, "    • Maximum relative deviation: $max_var_dev_rel")
  println(io, "    • Maximum relative deviation among observable variables: $max_observable_var_dev_rel")
  println(
    io,
    "    • Maximum relative deviation among variables with absolute deviation larger than $(solcomp.abstol): $max_var_dev_rel_ess",
  )
  println(
    io,
    "    • Maximum relative deviation among observable variables with absolute deviation larger than $(solcomp.abstol): $max_observable_var_dev_rel_ess",
  )
  println(io, "  • Equations...")
  if isequal(neqs1, neqs2)
    println(io, "    • ... in underlying problem: $neqs1")
  else
    println(io, "    • ... in underlying problem for the first solution: $neqs1")
    println(io, "    • ... in underlying problem for the second solution: $neqs2")
  end
  println(io, "    • ... shared by both solutions: $neqs_shared")
  println(io, "    • Maximum absolute deviation: $max_eq_dev_abs")
end

"""
Compares two provided ReconciliationSolutions.
It is assumed that the underlying ReconciliationProblems are identical. (This is roughly checked by comparing the problem names, the variable names and the types of shared variables.)
For each variable that has assigned values in both ReconciliationSolutions, the absolute difference of the assigned values and the relative difference (absolute difference divided by the mean of the assigned values) is returned.
For each equation that has assigned values in both ReconciliationSolutions, the absolute difference of the assigned values is returned. (Note that two equations are considered to be the same only if they have exactly the same form, including the order of variables etc.)
"""
function compare(rsol1::ReconciliationSolution, rsol2::ReconciliationSolution)
  problem_name = rsol1.problem.name
  @assert isequal(rsol1.problem.name, rsol2.problem.name)
  shared_variables = intersect(keys(rsol1.variables_value), keys(rsol2.variables_value))
  @assert [rsol1.problem.variables[var] for var in shared_variables] ==
          [rsol2.problem.variables[var] for var in shared_variables]
  vars_absolute_deviation = Dictionary{String, Float64}()
  vars_relative_deviation = Dictionary{String, Float64}()
  for var in shared_variables
    varval1 = rsol1.variables_value[var]
    varval2 = rsol2.variables_value[var]
    absdiff = abs(varval1 - varval2)
    set!(vars_absolute_deviation, var, absdiff)
    if varval1 == 0.0 && varval2 == 0.0
      set!(vars_relative_deviation, var, 0.0)
    else
      set!(vars_relative_deviation, var, 2.0 * absdiff / (abs(varval1 + varval2)))
    end
  end

  @assert issetequal(keys(rsol1.problem.variables), keys(rsol2.problem.variables))
  vars_equal_observability = Dictionary{Label, Bool}()
  for var in keys(rsol1.problem.variables)
    equal_observability =
      (haskey(rsol1.unobservable_variables, var) && haskey(rsol2.unobservable_variables, var)) ||
      (!haskey(rsol1.unobservable_variables, var) && !haskey(rsol2.unobservable_variables, var))
    set!(vars_equal_observability, var, equal_observability)
  end

  shared_equations_ids = intersect(keys(rsol1.score.equations_value), keys(rsol2.score.equations_value))
  shared_equations = Dictionary(shared_equations_ids, [rsol1.problem.equations[eid] for eid in shared_equations_ids])
  eqs_absolute_deviation = Dictionary{EID, Float64}()
  for eid in shared_equations_ids
    absdiff = abs(rsol1.score.equations_value[eid] - rsol2.score.equations_value[eid])
    set!(eqs_absolute_deviation, eid, absdiff)
  end

  SolutionComparison(
    (rsol1, rsol2),
    problem_name,
    Dictionary{String, Type{<:Variable}}(shared_variables, [rsol1.problem.variables[var] for var in shared_variables]),
    vars_absolute_deviation,
    vars_relative_deviation,
    vars_equal_observability,
    shared_equations,
    eqs_absolute_deviation,
  )
end

# ==================
# Error propagation
# ==================

"""
Assign each measured, unmeasured and fixed variables (in that order) to an integer index starting at 1.
"""
function get_variable_indices(problem::ReconciliationProblem)::Dictionary{Label, Int64}
  vars = Label[]
  append!(vars, keys(problem.measured))
  append!(vars, keys(problem.unmeasured))
  append!(vars, keys(problem.fixed))
  idxs = 1:length(problem.variables)
  Dictionary(vars, idxs)
end

"Convenience alias for sparse matrices in compressed sparse column format."
const SparseMatrix = SparseMatrixCSC{Float64, Int64}

"""
Compute the Jacobians `Jx` and `Jy` for measured and unmeasured variables respectively.

The dictionary `variable_index` maps variable names to indices in the full matrix `[Jx Jy]`.
For efficiency, this method assumes that the problem has no fixed variables.
If that is not the case you should apply `remove_fixed_variables!` as pre-processing step.
"""
function jacobians(
  problem::ReconciliationProblem,
  vals::Dictionary{Label, Float64},
)::@NamedTuple{
  Jx::SparseMatrix,
  Jy::SparseMatrix,
  variable_index::Dictionary{Label, Int64},
}
  if !isempty(problem.fixed)
    throw(ArgumentError("The problem has fixed variables."))
  end

  variable_index = get_variable_indices(problem)

  # Jacobian of constraints c(x, y) wrt x
  Ix = Int[]
  Jx = Int[]
  Vx = Float64[]
  # Jacobian of constraints c(x, y) wrt y
  Iy = Int[]
  Jy = Int[]
  Vy = Float64[]

  nx = length(problem.measured)
  ny = length(problem.unmeasured)
  nc = length(problem.equations)

  for (c, eq) in enumerate(problem.equations)
    for t in eq.linear_terms
      @assert haskey(variable_index, t.name) "Variable $(t.name) not found."
      ind = variable_index[t.name]
      if ind <= nx
        push!(Ix, c)
        push!(Jx, ind)
        push!(Vx, t.factor)
      elseif ind <= nx + ny
        push!(Iy, c)
        push!(Jy, ind - nx)
        push!(Vy, t.factor)
      end
    end

    for t in eq.bilinear_terms
      # For each bilinear term `t * var1 * var2`, compute the derivative respect
      # to each variable.
      # The combine argument (`+`) in the sparse matrix constructor ensures that
      # elements in repeated entries are added.

      # Derivative wrt first variable.
      ind = variable_index[t.name1]
      if ind <= nx
        push!(Ix, c)
        push!(Jx, ind)
        push!(Vx, t.factor * vals[t.name2])
      elseif ind <= nx + ny
        push!(Iy, c)
        push!(Jy, ind - nx)
        push!(Vy, t.factor * vals[t.name2])
      end

      # Derivative wrt second variable.
      ind = variable_index[t.name2]
      if ind <= nx
        push!(Ix, c)
        push!(Jx, ind)
        push!(Vx, t.factor * vals[t.name1])
      elseif ind <= nx + ny
        push!(Iy, c)
        push!(Jy, ind - nx)
        push!(Vy, t.factor * vals[t.name1])
      end
    end
  end

  Jx = sparse(Ix, Jx, Vx, nc, nx, +)
  Jy = sparse(Iy, Jy, Vy, nc, ny, +)
  (; Jx, Jy, variable_index)
end

"""
Compute Jacobians and Hessian for a reconciliation problem.

The theoretical derivation can be found in `resources/error_propagation.qmd`.
"""
function jacobians_and_hessian(sol::ReconciliationSolution)
  λ = sol.infos[INFO_DUAL]::Vector{Float64}
  problem = sol.infos[INFO_REDUCED_PROBLEM]
  variables_value = sol.variables_value
  jacobians_and_hessian(problem, variables_value, λ)
end
function jacobians_and_hessian(
  problem::ReconciliationProblem,
  variables_value::Dictionary{Label, Float64},
  λ::Vector{Float64},
)
  variable_index = get_variable_indices(problem)
  # Jacobian of constraints c(x, y) wrt x
  Ix = Int[]
  Jx = Int[]
  Vx = Float64[]
  # Jacobian of constraints c(x, y) wrt y
  Iy = Int[]
  Jy = Int[]
  Vy = Float64[]
  # Hessian of λ^T * c(x, y) wrt x and y
  IH = Int[]
  JH = Int[]
  VH = Float64[]
  nx = length(problem.measured)
  ny = length(problem.unmeasured)
  nc = length(problem.equations)
  for (c, eq) in enumerate(problem.equations)
    for t in eq.linear_terms
      ind = variable_index[t.name]
      if ind <= nx
        push!(Ix, c)
        push!(Jx, ind)
        push!(Vx, t.factor)
      elseif ind <= nx + ny
        push!(Iy, c)
        push!(Jy, ind - nx)
        push!(Vy, t.factor)
      end
    end
    for t in eq.bilinear_terms
      ind1 = variable_index[t.name1]
      ind2 = variable_index[t.name2]
      if ind1 <= nx
        push!(Ix, c)
        push!(Jx, ind1)
        push!(Vx, t.factor * variables_value[t.name2])
      elseif ind1 <= nx + ny
        push!(Iy, c)
        push!(Jy, ind1 - nx)
        push!(Vy, t.factor * variables_value[t.name2])
      end
      if ind2 <= nx
        push!(Ix, c)
        push!(Jx, ind2)
        push!(Vx, t.factor * variables_value[t.name1])
      elseif ind2 <= nx + ny
        push!(Iy, c)
        push!(Jy, ind2 - nx)
        push!(Vy, t.factor * variables_value[t.name1])
      end
      if ind1 <= nx + ny && ind2 <= nx + ny
        push!(IH, ind1)
        push!(JH, ind2)
        push!(IH, ind2)
        push!(JH, ind1)
        push!(VH, t.factor * λ[c])
        push!(VH, t.factor * λ[c])
      end
    end
  end
  return (;
    Jx = sparse(Ix, Jx, Vx, nc, nx),
    Jy = sparse(Iy, Jy, Vy, nc, ny),
    H = sparse(IH, JH, VH, nx + ny, nx + ny),
  )
end

function propagate_errors(sol::ReconciliationSolution; qr_tol = 1e-14, pinv_method::Symbol = :qr)
  (; Jx, Jy, H) = jacobians_and_hessian(sol)
  (; measured, unmeasured, equations) = sol.infos[INFO_REDUCED_PROBLEM]
  measured_names = values(map(x -> x.name, measured)) |> collect
  unmeasured_names = values(map(x -> x.name, values(unmeasured))) |> collect
  nx = length(measured)
  if nx == 0
    @debug "There are no measured variables. No errors to propagate."
    names = vcat(measured_names, unmeasured_names)
    return Dictionary{Label, Float64}()
  end
  ny = length(unmeasured)
  nc = length(equations)
  variances = zeros(nx)
  for (i, m) in enumerate(measured)
    variances[i] = m.uncertainty^2
  end
  Σ = Diagonal(variances)
  Hxx = H[1:nx, 1:nx]
  Hxy = H[1:nx, (nx + 1):(nx + ny)]
  Hyy = H[(nx + 1):(nx + ny), (nx + 1):(nx + ny)]
  K = [
    (2 * sparse(inv(Σ))+Hxx) Hxy Jx';
    Hxy' Hyy Jy';
    Jx Jy spzeros(nc, nc)
  ]
  fact = qr(K, tol = qr_tol)
  pinvK = zeros(size(K))
  b = zeros(size(K, 1))
  for j in 1:size(K, 2)
    if j > 1
      b[j - 1] = 0
    end
    b[j] = 1
    pinvK[:, j] = fact \ b
  end
  pinvK_x = pinvK[1:nx, 1:nx]
  iΣ = inv(Σ)
  x_std = diag(4 * pinvK_x * (iΣ * pinvK_x')) .|> sqrt
  if ny == 0
    y_std = Float64[]
  else
    pinvK_y = pinvK[(nx + 1):(nx + ny), 1:nx]
    y_std = diag(4 * pinvK_y * (iΣ * pinvK_y')) .|> sqrt
  end
  names = vcat(measured_names, unmeasured_names)
  stds = vcat(x_std, y_std)
  # If uncertainty information is not available it should be discarded.
  result = Dictionary(names, stds)
  filter!(!isnan, result)
end

function propagate_errors_RREF(sol::ReconciliationSolution)
  if haskey(sol.infos, INFO_REDUCED_PROBLEM)
    rec = sol.infos[INFO_REDUCED_PROBLEM]
  else
    rec = sol.problem
  end
  measured_vars = [k for k in keys(rec.measured)]
  unmeasured_vars = [k for k in keys(rec.unmeasured)]
  num_x = length(measured_vars)
  if num_x == 0
    @debug "There are no measured variables. No errors to propagate."
    return Dictionary{Label, Float64}()
  end
  num_y = length(unmeasured_vars)
  sol_vals = sol.variables_value
  variances = zeros(num_x)
  for (i, m) in enumerate(rec.measured)
    variances[i] = m.uncertainty^2
  end
  Q = Diagonal(variances)
  (; Jx, Jy, Jz, b) = jacobian(rec, sol_vals) # linearize problem at solution values
  Jyx = hcat(Jy, Jx)
  rref, pivots = rref_with_pivots!(Jyx, DEFAULT_ZTOL) # first the unmeasured then the measured variables to obtain equations only containing measured variables
  constr_after_y = findfirst(x -> x > num_y, pivots) # index of the first row with zero entries for unmeasured variables
  num_constr_y = (isnothing(constr_after_y) ? length(pivots) : constr_after_y - 1) # number of rows with nonzero entries for unmeasured variables
  num_constr_x = length(pivots) - num_constr_y # number of rows with zero entries for unmeasured variables and nonzero entries for measured variables
  # Identify unobservable variables and corresponding equations
  col_to_delete = Set{Int64}()
  row_to_delete = Set{Int64}()
  curr_pivot = 0
  for r in 1:num_constr_y
    prev_pivot = curr_pivot
    curr_pivot = pivots[r]
    if curr_pivot - prev_pivot > 1
      union!(col_to_delete, (prev_pivot + 1):(curr_pivot - 1))
    end
    row = rref[r, 1:num_y] # all y-entries for a row
    nonzero_entries = findall(!isapproxzero, row)
    if length(nonzero_entries) > 1
      push!(row_to_delete, r)
      for i in nonzero_entries
        push!(col_to_delete, i)
      end
    end
  end
  if !isequal(curr_pivot, num_y)
    union!(col_to_delete, (curr_pivot + 1):num_y)
  end
  # Remove rows and columns for unobservable variables
  col_to_keep = filter(i -> i ∉ col_to_delete, axes(rref, 2))
  row_to_keep = filter(i -> i ∉ row_to_delete, axes(rref, 1))
  essential_rref = rref[row_to_keep, col_to_keep]
  unmeasured_vars = unmeasured_vars[filter(i -> i <= num_y, col_to_keep)]
  num_y += -length(col_to_delete)
  num_constr_y += -length(row_to_delete)
  @assert num_y == num_constr_y
  # We decompose essential_rref according to the STAN paper
  Acy = essential_rref[1:num_constr_y, 1:num_y]
  Acx = essential_rref[1:num_constr_y, (num_y + 1):(num_y + num_x)]
  Arx = essential_rref[(num_constr_y + 1):(num_constr_y + num_constr_x), (num_y + 1):(num_y + num_x)]
  @assert isapprox(Acy, I) # after removing unobservable variables, Acy should be the identity matrix
  # Compute variance-covariance matrices
  temp_mat = (Arx * Q * Arx') \ (Arx * Q)
  cov_x = Q - (Q * Arx') * temp_mat
  cov_y = Acx * cov_x * Acx'
  uncertainties_x = diag(cov_x) .|> sqrt
  uncertainties_y = diag(cov_y) .|> sqrt
  Dictionary(vcat(measured_vars, unmeasured_vars), vcat(uncertainties_x, uncertainties_y))
end


"""
Computes Jacobian matrix for a ReconciliationProblem at a given point p as a dense matrix.
It returns matrices (; Jx, Jy, Jz, b), where each row is an equation and each column is a variable (measured in Jx, unmeasured in Jy, fixed in Jz) or the constant terms (b), 
i.e., the equation system is approximated by the linearization [Jx Jy Jz] * [x^T, y^T, z^T]^T = -b for measured variables x, unmeasured variables y and fixed variables z.
Computing the Jacobian matrix for bilinear terms requires the evaluation at some p. If p doesn't provide the evaluation values for a bilinear term, the equation including it is skipped and not part of Jx, Jy, Jz and b.
"""
function jacobian(rec::ReconciliationProblem, p::Dictionary{String, Float64} = Dictionary{String, Float64}())
  variable_index = get_variable_indices(rec)
  num_measured = length(rec.measured)
  num_unmeasured = length(rec.unmeasured)
  num_fixed = length(rec.fixed)
  num_equations = length(rec.equations)
  Jx = zeros(Float64, num_equations, num_measured)
  Jy = zeros(Float64, num_equations, num_unmeasured)
  Jz = zeros(Float64, num_equations, num_fixed)
  b = zeros(Float64, num_equations)
  i = 1 # counter for current row in jacobian
  warnings = Dictionary{EID, String}() # Dictionary to store warnings instead of showing all
  for eqid in keys(rec.equations)
    eq = rec.equations[eqid]
    valid_equation = true # indicates if the equation could be differentiated and evaluated successfully
    # Derivative of linear terms
    for t in eq.linear_terms
      var = t.name
      if haskey(rec.measured, var)
        j = variable_index[var]
        Jx[i, j] += t.factor
      elseif haskey(rec.unmeasured, var)
        j = variable_index[var] - num_measured
        Jy[i, j] += t.factor
      elseif haskey(rec.fixed, var)
        j = variable_index[var] - num_measured - num_unmeasured
        Jz[i, j] += t.factor
      else
        throw(ArgumentError("The variable $var does not exist in the problem."))
      end
    end
    # Derivative of bilinear terms
    for t in eq.bilinear_terms
      var1 = t.name1
      var2 = t.name2
      # The case that evaluation of the bilinear term is not possible due to missing information in p
      if !haskey(p, var1) || !haskey(p, var2)
        set!(
          warnings,
          eqid,
          "Warning for equation $eqid: The point at which linearization shall performed is not sufficiently defined (value for $(t.name1) or $(t.name2) is missing). Equation $eq deleted from matrix. This can result in a loss of information, potentially causing variables to be incorrectly identified as ambiguous.",
        )
        # Delete the current row
        Jx = Jx[1:end .!= i, :]
        Jy = Jy[1:end .!= i, :]
        Jz = Jz[1:end .!= i, :]
        b = b[1:end .!= i]
        valid_equation = false
        break # skip checking the remaining bilinear terms
      end
      # Derivatives of c*x*y with respect to the components x and y evaluated at p
      # Derivative wrt the first variable var1
      if haskey(rec.measured, var1)
        j1 = variable_index[var1]
        Jx[i, j1] += t.factor * p[var2]
      elseif haskey(rec.unmeasured, var1)
        j1 = variable_index[var1] - num_measured
        Jy[i, j1] += t.factor * p[var2]
      elseif haskey(rec.fixed, var1)
        j1 = variable_index[var1] - num_measured - num_unmeasured
        Jz[i, j1] += t.factor * p[var2]
      else
        throw(ArgumentError("The variable $var1 does not exist in the problem."))
      end
      # Derivative wrt the second variable var2
      if haskey(rec.measured, var2)
        j2 = variable_index[var2]
        Jx[i, j2] += t.factor * p[var1]
      elseif haskey(rec.unmeasured, var2)
        j2 = variable_index[var2] - num_measured
        Jy[i, j2] += t.factor * p[var1]
      elseif haskey(rec.fixed, var2)
        j2 = variable_index[var2] - num_measured - num_unmeasured
        Jz[i, j2] += t.factor * p[var1]
      else
        throw(ArgumentError("The variable $var2 does not exist in the problem."))
      end
    end
    # if the equation is not valid (i.e., the row was deleted), skip the rest of the for-loop
    if !valid_equation
      continue
    end
    # Store constant terms
    b[i] = eq.constant_term.value
    i += 1 # go to next row if equation was stored correctly
  end
  # Warn user if some equations were deleted
  if !iszero(length(warnings))
    @warn "$(length(warnings)) equations have been deleted from the Jacobian matrix. This can result in a loss of information, potentially causing variables to be incorrectly identified as ambiguous."
  end
  (; Jx, Jy, Jz, b)
end

end
