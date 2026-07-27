"""
Transformations that simplify a reconciliation problem before optimization.

The presolver applies a sequence of structure-exploiting steps:
- Trivial substitutions: solve for variables that appear in single-variable equations.
- Bilinear substitutions: propagate known transfer coefficients into adjacent equations.
- Fixed variable removal: substitute fixed values into all equations.
- Redundant equation removal: detect and drop linearly dependent equations using RREF or QR.
- Unobservable variable detection: identify unmeasured variables not uniquely determined by the system.

The result is a `PresolvedReconciliationProblem` containing a smaller equivalent reduced
problem plus a dictionary of already-determined variable assignments.
"""
module PresolveModule

using Reexport
using Accessors
using Dictionaries
using LinearAlgebra
using SparseArrays
using ResultTypes: Result, @try, iserror, unwrap

using ..ConstraintsModule
using ..ComparisonModule
using ..ErrorModule
using ..ReconciliationModule
using ..RowEchelonModule
using ..SparseRowEchelonModule

export presolve,
  presolve!,
  replace_variables!,
  replace_bilinear_terms!,
  remove_fixed_variables!,
  remove_measured_variables!,
  trivial_substitutions!,
  trivial_substitutions_recursive!,
  bilinear_substitutions!,
  trivial_and_bilinear_substitutions_recursive!,
  PresolvedReconciliationProblem,
  RedundancyData,
  RREFConfig,
  SparseRREFConfig,
  QRConfig,
  PresolveConfig,
  update_unobservable_vars!,
  remove_redundant_equations!,
  equations_matrix,
  independent_cols_qr,
  independent_cols!,
  augmented_matrix,
  find_unobservable_variables_in_constraints!,
  find_unobservable_variables_in_constraints

"For now only refers to the linear dependency solver."
abstract type PresolveConfig end

"Config using RREF for linear dependency reduction."
@kwdef struct RREFConfig <: PresolveConfig
  "Whether to scale the equations matrix before RREF."
  scale::Bool = true
  "Absolute tolerance for zero detection."
  atol::Float64 = DEFAULT_ZTOL
end

"Config using sparse RREF for linear dependency reduction."
@kwdef struct SparseRREFConfig <: PresolveConfig
  "Whether to scale the equations matrix before sparse RREF."
  scale::Bool = true
  "Absolute tolerance for zero detection."
  atol::Float64 = DEFAULT_ZTOL
end

"Config using QR for linear dependency reduction."
@kwdef struct QRConfig <: PresolveConfig
  "Whether to scale the equations matrix before QR decomposition."
  scale::Bool = true
  "Relative tolerance for rank detection. If nothing, it is computed based on the equations matrix."
  rtol::Union{Float64, Nothing} = nothing
  "Whether to use sparse QR decomposition. If nothing, it is computed based on the equations matrix."
  sparse::Union{Bool, Nothing} = nothing
end

"Use RREF for linear dependency reduction depending on sparsity pattern."
@kwdef struct MixConfig <: PresolveConfig
  "Whether to scale the equations matrix before reduction."
  scale::Bool = true
  "Configuration for dense RREF (used when sparsity is below threshold)."
  rref_config::RREFConfig = RREFConfig()
  "Configuration for sparse RREF (used when sparsity is above threshold)."
  sparse_config::SparseRREFConfig = SparseRREFConfig()
  "Sparsity threshold (fraction of zeros) to choose between dense and sparse RREF."
  sparsity::Float64 = 0.8
end

"Default presolve config."
PresolveConfig() = MixConfig()


"""
Data structure holding the result of redundant equations removal.
"""
struct RedundancyData
  "Pivots matrix from the rref simplification routine."
  pivots::Vector{Int64}
  "Constant terms vector (one entry per equation)."
  c::Vector{Float64}
  "Linear/bilinear terms matrix (rows = terms, columns = equations)."
  A::Matrix{Float64}
  "Indices (identifier) of equations which were removed."
  removed_eqs_eid::Vector{EID}
  "Indices of linear terms in the equations matrix."
  linear_term_idx::Dictionary{Label, Int64}
  "Indices of bilinear terms in the equations matrix."
  bilinear_term_idx::Dictionary{Tuple{Label, Label}, Int64}
end

"""
Struct holding the problem after the symbolic/numeric presolve pipeline.
"""
@kwdef mutable struct PresolvedReconciliationProblem
  "Resulting reconciliation problem."
  pre::ReconciliationProblem
  "Variables which have been replaced and their corresponding values."
  replacements::Dictionary{Label, Float64} = Dictionary{Label, Float64}()
  "Variable names which didn't get replaced. These are remaining variable to solve for after presolve."
  variables_kept::Vector{Label}
  "Variable names which got replaced by specific values."
  variables_removed::Vector{Label} = Vector{Label}()
  "Reduced redundant equations data structure."
  rdata::Union{RedundancyData, Nothing} = nothing
  "Measured variables that don't occur in constraint equations or replacements after presolving."
  unconstrained::Dictionary{Label, MeasuredVariable} = Dictionary{Label, MeasuredVariable}()
  "Unmeasured variables that don't occur in constraint equations or replacements after presolving."
  unobservable::Dictionary{Label, UnmeasuredVariable} = Dictionary{Label, UnmeasuredVariable}()
end

function PresolvedReconciliationProblem(
  rec::ReconciliationProblem;
  replacements::Dictionary{Label, Float64} = Dictionary{Label, Float64}(),
  variables_removed::Vector{Label} = Vector{Label}(),
  rdata::Union{RedundancyData, Nothing} = nothing,
  unconstrained::Dictionary{Label, MeasuredVariable} = Dictionary{Label, MeasuredVariable}(),
  unobservable::Dictionary{Label, UnmeasuredVariable} = Dictionary{Label, UnmeasuredVariable}(),
)
  PresolvedReconciliationProblem(;
    pre = rec,
    variables_kept = collect(keys(rec.variables)),
    replacements,
    variables_removed,
    rdata,
    unconstrained,
    unobservable,
  )
end

"""
Check for variables in 'vars_to_check' if they don't participate in any equation of the PresolvedReconciliationProblem 'prerec'.
If so, declare it as 'unobservable' and remove it from 'rec'.
If 'vars_to_check' is not specified (i.e., 'nothing'), check all variables from the ReconciliationProblem.
"""
function update_unobservable_vars!(
  prerec::PresolvedReconciliationProblem,
  vars_to_check::Union{Vector{Label}, Nothing} = nothing,
)
  rec = prerec.pre
  if isnothing(vars_to_check)
    vars_to_check = collect(keys(rec.variables))
  end
  # Remove variables from the problem if there are no equations in which they participate.
  for name in vars_to_check
    if !haskey(rec.variables, name)
      # Skip if already removed.
      continue
    end
    # If there are no more equations in which this variable participates, then remove the variable.
    idx = rec.var_to_eq[name]
    if isempty(idx) && !haskey(rec.fixed, name)
      if haskey(rec.measured, name)
        insert!(prerec.unconstrained, name, rec.measured[name])
      elseif haskey(rec.unmeasured, name)
        insert!(prerec.unobservable, name, rec.unmeasured[name])
      end
      push!(prerec.variables_removed, name)
      filter!(var -> var != name, prerec.variables_kept)
      remove_variable!(rec, name)
      @debug "Variable $name has been removed as a side effect of replacements."
    end
  end
  prerec
end
function update_unobservable_vars!(rec::ReconciliationProblem, vars_to_check::Union{Vector{Label}, Nothing} = nothing)
  prerec = PresolvedReconciliationProblem(rec)
  update_unobservable_vars!(prerec, vars_to_check)
end

"""
Each given variable is replaced by its corresponding constant value and removed from the problem.
The result is an equivalent problem without those variables.

Important: The `variables` dictionary cannot be alised to any of `rec`'s dictionaries.

## Algorithm

For each variable we loop through each equation, making the replacements:

- Linear terms with matching name are replaced with constants.
- Bilinear terms with matching name are replaced to linear terms.

For the tolerance for the consistency check of equations (left-hand side is zero) during variable removal
see `ComparisonModule`.
"""
function replace_variables!(
  prerec::PresolvedReconciliationProblem,
  variables::Dictionary{Label, Float64};
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  rec = prerec.pre
  # Variables `x` appearing in terms of the form `x * v` where v has been replaced by 0
  # can potentially be removed from the problem as there might not be other equations in which `x` participate.
  maybe_remove_from_problem = Label[]

  # Storing replacements done for error logging
  replacements = Dictionary{Label, Float64}()

  # Replace the variable on each term in which it participates.
  # Linear terms turn into constant terms, while bilinear terms turn into linear terms.
  for (n, v) in pairs(variables)
    set!(replacements, n, v)

    # Remove variable from unconstrained and unobservable variables
    if haskey(prerec.unconstrained, n)
      delete!(prerec.unconstrained, n)
    end
    if haskey(prerec.unobservable, n)
      delete!(prerec.unobservable, n)
    end

    # Only iterate over equations in which this variable participates.
    for eid in rec.var_to_eq[n]
      maybe_remove_from_equation = Label[]
      eq = rec.equations[eid]
      for lt in eq.linear_terms
        if lt.name == n
          val = lt.factor * v
          if !isapproxzero(val)
            ct = ConstantTerm(val)
            eq.constant_term += ct
          end
        end
      end
      # Only keep linear terms which differ from the replacement.
      filter!(lt -> lt.name != n, eq.linear_terms)

      for bt in eq.bilinear_terms
        is_name1 = n == bt.name1
        is_name2 = n == bt.name2
        if is_name1 || is_name2
          if is_name1 && is_name2
            # The matching term is of the form `var^2`, so we push a constant.
            val = bt.factor * v^2
            eq.constant_term += ConstantTerm(val)
          else
            val = bt.factor * v
            # We replace the found variable, keeping the other one.
            other_name = is_name1 ? bt.name2 : bt.name1
            if isapproxzero(val)
              if other_name ∉ maybe_remove_from_equation
                push!(maybe_remove_from_equation, other_name)
              end
              continue
            end

            add_linear_term!(rec, eid, other_name, val, maybe_remove_from_equation)
          end
        end
      end
      # Only keep bilinear terms which differ from the replacement.
      filter!(bt -> (bt.name1 != n) && (bt.name2 != n), eq.bilinear_terms)
      # This variable does no longer participate in equation `eid` so we remove it from the reverse lookup dictionary.
      deleteat!(rec.eq_to_var[eid], findfirst(==(n), rec.eq_to_var[eid]))
      # Remove variables from var_to_eq and eq_to_var if they don't participate anymore
      update_vars_dicts!(rec, maybe_remove_from_equation, [eid], maybe_remove_from_problem)

      # Remove any spurious equation which might have been generated after replacements.
      # If an inconsistency is found, an error is produced.
      # TODO Improve by returning equations which were inconsistent down the pipeline.
      if isempty(eq.linear_terms) && isempty(eq.bilinear_terms)
        # By repeated variable replacements we may find an equation with only constants.
        # In that case check that the constant is zero.
        s = evaluate(eq.constant_term)
        if !(abs(s) < tol)
          return InfeasibleEquationError(
            "Equation $eid : $eq expected to be zero, evaluated to: $s.",
            [eid],
            [s],
            Dictionary(
              ["Reason for infeasibility", "Problematic replacements"],
              ["Contradicting equations after replacements", replacements],
            ),
          )
        end
        remove_equation!(rec, eid)
      end
    end

    # Since this variable has been replaced on each equation, it can be safely removed.
    push!(prerec.variables_removed, n)
    filter!(var -> var != n, prerec.variables_kept)
    remove_variable!(rec, n)
  end

  # Remove variables from the problem if there are no equations in which they participate.
  update_unobservable_vars!(prerec, maybe_remove_from_problem)

  merge!(prerec.replacements, replacements)
  prerec
end

function replace_variables!(
  rec::ReconciliationProblem,
  variables::Dictionary{Label, Float64};
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  prerec = PresolvedReconciliationProblem(rec)
  replace_variables!(prerec, variables; tol)
end


"""
Each given product `v * w` is replaced by its corresponding constant value and removed from the problem.
The result is an equivalent problem without those product terms provided that one equation defining the product value is kept (to be specified in storing_equations).
The method does not check for `w * v`!

Important: The `variables` dictionary cannot be alised to any of `rec`'s dictionaries.

## Algorithm

For each variable we loop through each equation, making the replacements:

- Bilinear terms with matching product (permutations are not considered!) are replaced to constant terms.

For the tolerance for the consistency check of equations (left-hand side is zero) during variable removal
see `ComparisonModule`.
"""
function replace_bilinear_terms!(
  prerec::PresolvedReconciliationProblem,
  terms::Dictionary{Tuple{Label, Label}, Float64},
  storing_equations::Dictionary{Tuple{Label, Label}, EID} = Dictionary{Tuple{Label, Label}, EID}();
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  rec = prerec.pre

  # Storing replacements done for error logging
  replacements = Dictionary{Tuple{Label, Label}, Float64}()

  # Replace the product on each term in which it participates.
  # The bilinear terms turn into constant terms.
  for ((n1, n2), v) in pairs(terms)
    set!(replacements, (n1, n2), v)
    # Check if the product with permutated factors is defined in the same way
    if !haskey(terms, (n2, n1)) || !relaxed_isapprox(terms[(n2, n1)], v)
      @warn "The products $n1 * $n2 and $n2 * $n1 are not defined in the same way."
    end
    # Only iterate over equations in which both variables participates.
    equation_list = intersect(rec.var_to_eq[n1], rec.var_to_eq[n2])
    for eid in equation_list
      # Skip equation if it is the one storing the product value
      if haskey(storing_equations, (n1, n2)) && eid == storing_equations[(n1, n2)]
        continue
      end
      eq = rec.equations[eid]

      for bt in eq.bilinear_terms
        if bt.name1 == n1 && bt.name2 == n2
          val = bt.factor * v
          # We replace the found product.
          if !iszero(val)
            ct = ConstantTerm(val)
            eq.constant_term += ct
          end
        end
      end
      # Only keep bilinear terms which differ from the replacement.
      filter!(bt -> (bt.name1, bt.name2) != (n1, n2), eq.bilinear_terms)
      # Check if the variables occur in the terms. If not, update the lookup dictionaries.
      for var in [n1, n2]
        var_found = false
        for lt in eq.linear_terms
          if lt.name == var
            var_found = true
          end
        end
        for bt in eq.bilinear_terms
          if bt.name1 == var || bt.name2 == var
            var_found = true
          end
        end
        if !var_found
          # Update lookup dictionaries
          deleteat!(rec.eq_to_var[eid], findall(x -> x == var, rec.eq_to_var[eid]))
          deleteat!(rec.var_to_eq[var], findall(x -> x == eid, rec.var_to_eq[var]))
        end
      end

      # Remove any spurious equation which might have been generated after replacements.
      # If an inconsistency is found, an error is produced.
      # TODO Improve by returning equations which were inconsistent down the pipeline.
      if isempty(eq.linear_terms) && isempty(eq.bilinear_terms)
        # By repeated variable replacements we may find an equation with only constants.
        # In that case check that the constant is zero.
        s = evaluate(eq.constant_term)
        if !(abs(s) < tol)
          return InfeasibleEquationError(
            "Equation $eid : $eq expected to be zero, evaluated to: $s.",
            [eid],
            [s],
            Dictionary(
              ["Reason for infeasibility", "Problematic product replacements"],
              ["Contradicting equations after product replacements", replacements],
            ),
          )

        end
        remove_equation!(rec, eid)
      end
    end
  end

  prerec
end

function replace_bilinear_terms!(
  rec::ReconciliationProblem,
  terms::Dictionary{Tuple{Label, Label}, Float64},
  storing_equations::Dictionary{Tuple{Label, Label}, EID} = Dictionary{Tuple{Label, Label}, EID}();
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  prerec = PresolvedReconciliationProblem(rec)
  replace_bilinear_terms!(prerec, terms, storing_equations; tol)
end


"""
Remove fixed variables replacing by their corresponding values.
"""
function remove_fixed_variables!(
  prerec::PresolvedReconciliationProblem;
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  fixed_vals = dictionary(k => v.value for (k, v) in pairs(prerec.pre.fixed))
  result = replace_variables!(prerec, fixed_vals; tol)
  if iserror(result)
    add_info!(result.error.infos, "Failure point", "Removal of fixed variables")
  else
    @assert isempty(prerec.pre.fixed)
  end
  return result
end
function remove_fixed_variables!(
  rec::ReconciliationProblem;
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  prerec = PresolvedReconciliationProblem(rec)
  remove_fixed_variables!(prerec; tol)
end

"""
Each measured variable is replaced by its corresponding mean value and removed from the problem.
The result is an equivalent problem without measured variables.
"""
function remove_measured_variables!(
  prerec::PresolvedReconciliationProblem;
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  measured_vals = dictionary(k => v.value for (k, v) in pairs(prerec.pre.measured))
  unconstrained_vals = dictionary(k => v.value for (k, v) in pairs(prerec.unconstrained))
  merge!(measured_vals, unconstrained_vals)
  result = replace_variables!(prerec, measured_vals; tol)
  if iserror(result)
    add_info!(result.error.infos, "Failure point", "Removal of measured variables")
  else
    @assert isempty(prerec.pre.measured)
    @assert isempty(prerec.unconstrained)
  end
  return result
end
function remove_measured_variables!(
  rec::ReconciliationProblem;
  tol::Float64 = DEFAULT_ZTOL,
)::Result{PresolvedReconciliationProblem, InfeasibleEquationError}
  prerec = PresolvedReconciliationProblem(rec)
  remove_measured_variables!(prerec; tol)
end

"""
Check if new replacements from trivial replacements are valid.
"""
function check_replacements!(
  rec::ReconciliationProblem,
  eid::EID,
  eq::Equation,
  replacement_keys::Vector{Label},
  val::Float64,
  rtol::Float64,
  new_replacements::Dictionary{Label, Float64},
  replacements::Dictionary{Label, Float64},
  replacement_equations::Dictionary{Label, EID},
  invalid_replacements::Dictionary{Label, Vector{Float64}},
  invalid_equations::Dictionary{Label, EID},
)
  # We cannot use `insert!` here since a replacement is not guaranteed to be unique:
  # for example, when replacing measured variables by their mean values it can happen that the same unmeasured variable
  # can be solved in two contradictory ways in two different equations.
  if haskey(new_replacements, replacement_keys[1]) || haskey(replacements, replacement_keys[1])
    previous_replacements = haskey(new_replacements, replacement_keys[1]) ? new_replacements : replacements
    prev = previous_replacements[replacement_keys[1]]
    if !relaxed_isapprox(val, prev; rtol)
      # If previous_replacements == replacements, this case should not occur since after (successful) previous replacements, only one equation defining t.name1 * t.name2 is left.
      # Value mismatch: invalid replacement and prevent adding the variable again.
      for repl_key in replacement_keys
        delete!(previous_replacements, repl_key)
        delete!(replacement_equations, repl_key)
        set!(invalid_replacements, repl_key, [prev, val])
        set!(invalid_equations, repl_key, eid)
      end
    end
  elseif haskey(invalid_replacements, replacement_keys[1])
    for repl_key in replacement_keys
      push!(invalid_replacements[repl_key], val)
    end
  else
    for repl_key in replacement_keys
      var = get_variable(rec, repl_key)
      if !isnothing(var) && !satisfies_bounds(var, val)
        set!(invalid_replacements, repl_key, [val])
        return InfeasibleEquationError(
          "Equation $eid : $eq uniquely defines variable $(repl_key) to be value $val but this violates the prescribed bounds.",
          [eid],
          [0.0],
          Dictionary(
            ["Reason for infeasibility", "Conflicting replacements", "Previous replacements"],
            ["Replacement conflicts variable bounds", invalid_replacements, replacements]),
        )
      else
        set!(new_replacements, repl_key, val)
        set!(replacement_equations, repl_key, eid)
      end
    end
  end
end

"""
Check if new replacements from bilinear replacements are valid.
"""
function check_replacements!(
  rec::ReconciliationProblem,
  eid::EID,
  eq::Equation,
  replacement_keys::Vector{Tuple{Label, Label}},
  val::Float64,
  rtol::Float64,
  new_replacements::Dictionary{Tuple{Label, Label}, Float64},
  replacements::Dictionary{Tuple{Label, Label}, Float64},
  replacement_equations::Dictionary{Tuple{Label, Label}, EID},
  invalid_replacements::Dictionary{Tuple{Label, Label}, Vector{Float64}},
  invalid_equations::Dictionary{Tuple{Label, Label}, EID},
)
  # We cannot use `insert!` here since a replacement is not guaranteed to be unique:
  # for example, when replacing measured variables by their mean values it can happen that the same unmeasured variable
  # can be solved in two contradictory ways in two different equations.
  if haskey(new_replacements, replacement_keys[1]) || haskey(replacements, replacement_keys[1])
    previous_replacements = haskey(new_replacements, replacement_keys[1]) ? new_replacements : replacements
    prev = previous_replacements[replacement_keys[1]]
    if !relaxed_isapprox(val, prev; rtol)
      # If previous_replacements == replacements, this case should not occur since after (successful) previous replacements, only one equation defining t.name1 * t.name2 is left.
      # Value mismatch: invalid replacement and prevent adding the variable again.
      for repl_key in replacement_keys
        delete!(previous_replacements, repl_key)
        delete!(replacement_equations, repl_key)
        set!(invalid_replacements, repl_key, [prev, val])
        set!(invalid_equations, repl_key, eid)
      end
    end
  elseif haskey(invalid_replacements, replacement_keys[1])
    for repl_key in replacement_keys
      push!(invalid_replacements[repl_key], val)
    end
  else
    for repl_key in replacement_keys
      set!(new_replacements, repl_key, val)
      set!(replacement_equations, repl_key, eid)
    end
  end
end

function linear_root(factor::Float64, constant::Float64)::Float64
  # Normalize near-zero constant terms before division to avoid spurious
  # non-zero substitution values from floating-point rounding.
  c_value = isapproxzero(constant) ? zero(Float64) : constant
  val = -c_value / factor
  if isapproxzero(val)
    # To mitigate propagation of numerical errors we reset zero values.
    val = zero(Float64)
  end
  return val
end

"""
A trivial substitution is applied on terms of the form `c - v = 0` or `c + v² = 0` where `c` is a constant and `v` is a variable.

For squared terms, the substitution is only done if the variable bounds uniquely determine the solution,
i.e. one of `sqrt(c)` and `-sqrt(c)` is not within the variable bounds or they match within `rol` relative tolerance.

If we can solve for `v` in two or more equations, this is considered a valid substition only if the found values
match within `rol` relative tolerance.
"""
function trivial_substitutions!(
  prerec::PresolvedReconciliationProblem;
  rtol::Float64 = 1e-6,
)::Result{Bool, SolverError}
  rec = prerec.pre
  # Ensure that all linear equations that have a single linear term are simplified.
  # Assumptions:
  # - Equations of the form `c + x = 0`, ie. single linear term and constant.
  #   Note that multiple linear terms with the same label such as `x + 2x = 0` are not handled.
  # - If two equations are contradictory, eg. `c1 + x = 0` and `c2 + x = 0` with `c1 != c2`,
  #   then such replacement is skipped.
  new_replacements = Dictionary{Label, Float64}()
  to_remove_eid = Dictionary{Label, EID}()
  invalid_replacements = Dictionary{Label, Vector{Float64}}()
  invalid_equations = Dictionary{Label, EID}()
  for (eid, eq) in pairs(rec.equations)
    if isempty(eq.bilinear_terms) && isone(length(eq.linear_terms))
      t = only(eq.linear_terms)
      c = eq.constant_term
      if isapproxzero(t.factor)
        if !isapproxzero(c.value)
          return InfeasibleEquationError(
            "Equation $eid : $eq expected to be zero, evaluated to: $(c.value).",
            [eid],
            [c.value],
            Dictionary(
              ["Reason for infeasibility", "Failure point", "Previous replacements"],
              ["Contradicting equations", "Trivial substitutions", merge(prerec.replacements, replacements)],
            ),
          )
        end
        continue
      end
      val = linear_root(t.factor, c.value)
      if haskey(prerec.replacements, t.name)
        return InternalError(
          "The variable $(t.name) should not be present in the ReconciliationProblem anymore as it has already been replaced.",
        )
      end
      result = check_replacements!(
        rec,
        eid,
        eq,
        [t.name],
        val,
        rtol,
        new_replacements,
        prerec.replacements,
        to_remove_eid,
        invalid_replacements,
        invalid_equations,
      )
      if iserror(result)
        add_info!(result.infos, "Failure point", "Trivial substitutions")
        return result
      end
    elseif isempty(eq.linear_terms) && isone(length(eq.bilinear_terms))
      # Check for equations having a single bilinear term of the form c + x^2 = 0.
      t = only(eq.bilinear_terms)
      if isequal(t.name1, t.name2)
        c = eq.constant_term
        if isapproxzero(t.factor)
          if !isapproxzero(c.value)
            return InfeasibleEquationError(
              "Equation $eid : $eq expected to be zero, evaluated to: $(c.value).",
              [eid],
              [c.value],
              Dictionary(
                ["Reason for infeasibility", "Failure point", "Previous replacements"],
                ["Contradicting equations", "Trivial substitutions", merge(prerec.replacements, replacements)],
              ),
            )
          end
          continue
        end
        val = linear_root(t.factor, c.value)
        # Set val to square root
        if val < 0.0
          return InfeasibleEquationError(
            "Equation $eid : $eq claims that the square $(t.name1)^2 is the negative value $val.",
            [eid],
            [val],
            Dictionary(
              ["Reason for infeasibility", "Failure point", "Conflicting replacements", "Previous replacements"],
              [
                "Negative square",
                "Trivial substitutions",
                invalid_replacements,
                prerec.replacements,
              ]),
          )
        else
          val = sqrt(val)
        end
        # Check if negative square root would also be valid
        if satisfies_bounds(get_variable(rec, t.name1), -val)
          if !satisfies_bounds(get_variable(rec, t.name1), val)
            # If positive square root is not valid, use the negative square root for the replacement.
            val = -val
          elseif !relaxed_isapprox(val, -val; rtol)
            # If the positive and the negative square root are valid, don't do any replacement as it is not uniquely determined.
            continue
          end
        end
        if haskey(prerec.replacements, t.name1)
          return InternalError(
            "The variable $(t.name1) should not be present in the ReconciliationProblem anymore as it has already been replaced.",
          )
        end
        result = check_replacements!(
          rec,
          eid,
          eq,
          [t.name1],
          val,
          rtol,
          new_replacements,
          prerec.replacements,
          to_remove_eid,
          invalid_replacements,
          invalid_equations,
        )
        if iserror(result)
          add_info!(result.infos, "Failure point", "Trivial substitutions")
          return result
        end
      end
    end
  end
  if !isempty(invalid_replacements)
    prob_var = keys(invalid_replacements)
    prob_var_str = "$(join(prob_var, ", "))"
    index = [invalid_equations[var] for var in prob_var]
    residual = [invalid_replacements[var][1] - invalid_replacements[var][2] for var in prob_var]
    return InfeasibleEquationError(
      "Contradictions when replacing following variable(s): $prob_var_str",
      index,
      residual,
      Dictionary(
        ["Reason for infeasibility", "Failure point", "Conflicting replacements", "Previous replacements"],
        [
          "Contradicting equations after replacements",
          "Trivial substitutions",
          invalid_replacements,
          prerec.replacements,
        ],
      ),
    )
  end
  found = !isempty(new_replacements)

  if found
    # Remove from the equations the ones for which we performed substitutions.
    foreach(values(to_remove_eid)) do eid
      remove_equation!(rec, eid)
    end

    # Replace each found variable to the assigned value and remove the variable from the model.
    result = replace_variables!(prerec, new_replacements)
    if iserror(result)
      add_info!(result.error.infos, "Failure point", "Trivial substitutions")
      add_info!(result.error.infos, "Previous replacements", prerec.replacements)
      @try result
    end
  end

  found
end

function trivial_substitutions!(
  rec::ReconciliationProblem;
  rtol::Float64 = 1e-6,
)::Result{Bool, SolverError}
  prerec = PresolvedReconciliationProblem(rec)
  trivial_substitutions!(prerec; rtol)
end


"""
A bilinear substitution is applied on terms of the form `c - v * w = 0` where `c` is a constant and `v` and `w` are a variables.
In other equations, `v * w` or `w * v` is replaced by `c`.

If we can solve for `v * w` or `w * v` in two or more equations, this is considered a valid substition only if the found values
match within `rtol` relative tolerance.
"""
function bilinear_substitutions!(
  prerec::PresolvedReconciliationProblem,
  replacements::Dictionary{Tuple{Label, Label}, Float64} = Dictionary{Tuple{Label, Label}, Float64}();
  rtol::Float64 = 1e-6,
)::Result{Bool, SolverError}
  rec = prerec.pre

  new_replacements = Dictionary{Tuple{Label, Label}, Float64}()
  replacement_storing_eid = Dictionary{Tuple{Label, Label}, EID}() # The replacement of a product cannot be done in all equations but the initial equation `v * w = c` has to be kept.
  invalid_replacements = Dictionary{Tuple{Label, Label}, Vector{Float64}}()
  invalid_equations = Dictionary{Tuple{Label, Label}, EID}()
  for (eid, eq) in pairs(rec.equations)
    if isempty(eq.linear_terms) && isone(length(eq.bilinear_terms))
      t = only(eq.bilinear_terms)
      c = eq.constant_term
      if isapproxzero(t.factor)
        if !isapproxzero(t.value)
          return InfeasibleEquationError(
            "Equation $eid : $eq expected to be zero, evaluated to: $(c.value).",
            [eid],
            [c.value],
            Dictionary(
              ["Reason for infeasibility", "Failure point", "Previous product replacements", "Previous replacements"],
              [
                "Contradicting equations",
                "Bilinear substitutions",
                merge(replacements, new_replacements),
                prerec.replacements,
              ],
            ),
          )
        end
        continue
      end
      val = linear_root(t.factor, c.value)
      result = check_replacements!(
        rec,
        eid,
        eq,
        [(t.name1, t.name2), (t.name2, t.name1)],
        val,
        rtol,
        new_replacements,
        replacements,
        replacement_storing_eid,
        invalid_replacements,
        invalid_equations,
      )
      if iserror(result)
        add_info!(result.infos, "Failure point", "Bilinear substitutions")
        return result
      end
    end
  end
  if !isempty(invalid_replacements)
    prob_varpairs = keys(invalid_replacements)
    prob_varpairs_str = "$(join(prob_varpairs, ", "))"
    index = [invalid_equations[varpair] for varpair in prob_varpairs]
    residual = [invalid_replacements[varpair][1] - invalid_replacements[varpair][2] for varpair in prob_varpairs]
    return InfeasibleEquationError(
      "Contradictions when replacing following variable(s): $prob_varpairs_str",
      index,
      residual,
      Dictionary(
        [
          "Reason for infeasibility",
          "Failure point",
          "Conflicting replacements",
          "Previous product replacements",
          "Previous replacements",
        ],
        [
          "Contradicting equations after replacements",
          "Bilinear substitutions",
          invalid_replacements,
          merge(replacements, new_replacements),
          prerec.replacements,
        ],
      ),
    )
  end
  found = !isempty(new_replacements)

  if found
    # Merge changes into the output dictionary.
    merge!(replacements, new_replacements)

    # Replace each found variable to the assigned value and remove the variable from the model.
    result = replace_bilinear_terms!(prerec, new_replacements, replacement_storing_eid)

    if iserror(result)
      add_info!(result.error.infos, "Failure point", "Bilinear substitutions")
      add_info!(result.error.infos, "Previous product replacements", replacements)
      add_info!(result.error.infos, "Previous replacements", prerec.replacements)
    end
    @try result
  end

  found
end

function bilinear_substitutions!(
  rec::ReconciliationProblem,
  replacements::Dictionary{Tuple{Label, Label}, Float64} = Dictionary{Tuple{Label, Label}, Float64}();
  rtol::Float64 = 1e-6,
)::Result{Bool, SolverError}
  prerec = PresolvedReconciliationProblem(rec)
  bilinear_substitutions!(prerec, replacements; rtol)
end

"Recursively perform trivial substitions until no more are found."
function trivial_substitutions_recursive!(
  prerec::PresolvedReconciliationProblem;
  rtol::Float64 = 1e-6,
)::Result{PresolvedReconciliationProblem, SolverError}

  found = true
  while found
    # Iterate while there are replacements found.
    # If there's an infeasible equation, we can't continue with trivial substitutions.
    found = @try trivial_substitutions!(prerec; rtol)
  end

  prerec
end

function trivial_substitutions_recursive!(
  rec::ReconciliationProblem;
  rtol::Float64 = 1e-6,
)::Result{PresolvedReconciliationProblem, SolverError}
  prerec = PresolvedReconciliationProblem(rec)
  trivial_substitutions_recursive!(prerec; rtol)
end

"Recursively perform trivial and product substitions until no more are found."
function trivial_and_bilinear_substitutions_recursive!(
  prerec::PresolvedReconciliationProblem,
  product_replacements::Dictionary{Tuple{Label, Label}, Float64} = Dictionary{Tuple{Label, Label}, Float64}();
  rtol::Float64 = 1e-6,
)::Result{PresolvedReconciliationProblem, SolverError}

  found_bilinear = true
  while found_bilinear
    found_trivial = true
    while found_trivial
      # Iterate while there are trivial replacements found.
      # If there's an infeasible equation, we can't continue with trivial substitutions.
      trivial_result = trivial_substitutions!(prerec; rtol)
      if iserror(trivial_result) && !isempty(product_replacements)
        add_info!(trivial_result.error.infos, "Previous product replacements", product_replacements)
      end
      found_trivial = @try trivial_result

    end

    found_bilinear = @try bilinear_substitutions!(prerec, product_replacements; rtol)
  end

  prerec
end

function trivial_and_bilinear_substitutions_recursive!(
  rec::ReconciliationProblem,
  product_replacements::Dictionary{Tuple{Label, Label}, Float64} = Dictionary{Tuple{Label, Label}, Float64}();
  rtol::Float64 = 1e-6,
)::Result{PresolvedReconciliationProblem, SolverError}
  prerec = PresolvedReconciliationProblem(rec)
  trivial_and_bilinear_substitutions_recursive!(prerec, product_replacements; rtol)
end

"Scale in-place the columns of the matrix by the maximum absolute among rows. Returns the vector of scale factors applied to each column."
function scale_columns!(M)
  scales = ones(Float64, size(M, 2))
  for j in 1:size(M, 2)
    max_val = maximum(abs, @view(M[:, j]))
    if max_val > 0.0
      @. M[:, j] /= max_val
      scales[j] = max_val
    end
  end
  return M, scales
end

"""
Create matrix `A` (rows = variables, columns = equations) and constant terms vector `c`.
Bilinear terms are lifted to a new independent variable.
`linear_term_idx` maps variable labels to 1-based row indices in `A`.
`bilinear_term_idx` maps variable-pair labels to row indices in `A` continuing after the linear rows.
"""
function equations_matrix(rec::ReconciliationProblem)
  linear_term_idx = Dictionary{Label, Int64}()
  N = 0
  for group in (rec.measured, rec.unmeasured, rec.fixed)
    for (i, n) in enumerate(keys(group))
      N += 1
      set!(linear_term_idx, n, N)
    end
  end

  bilinear_term_idx = Dictionary{Tuple{Label, Label}, Int64}()
  N = length(linear_term_idx)
  for eq in rec.equations
    for t in eq.bilinear_terms
      if !haskey(bilinear_term_idx, (t.name1, t.name2)) && !haskey(bilinear_term_idx, (t.name2, t.name1))
        N += 1
        set!(bilinear_term_idx, (t.name1, t.name2), N)
      end
    end
  end

  neqs = length(rec.equations)
  nrows = length(linear_term_idx) + length(bilinear_term_idx)
  c = zeros(Float64, neqs)
  A = zeros(Float64, nrows, neqs)
  for (i, eq) in enumerate(rec.equations)
    c[i] = eq.constant_term.value
    for t in eq.linear_terms
      j = linear_term_idx[t.name]
      A[j, i] = t.factor
    end
    for t in eq.bilinear_terms
      j = if haskey(bilinear_term_idx, (t.name1, t.name2))
        bilinear_term_idx[(t.name1, t.name2)]
      else
        bilinear_term_idx[(t.name2, t.name1)]
      end
      A[j, i] = t.factor
    end
  end
  (; c, A, linear_term_idx, bilinear_term_idx)
end

"Return vector of pivots and nullspace of a matrix"
function independent_cols!(M::AbstractMatrix, config::RREFConfig)::Tuple{Vector{Int}, AbstractMatrix}
  rref, pivots = rref_with_pivots!(M, config.atol)

  ncols = size(rref, 2)
  nullvectors = Vector{Vector{Float64}}()
  for i in 1:ncols
    if i ∉ pivots
      # Column is dependent
      col = Vector(rref[:, i])
      col .= ifelse.(isapproxzero.(col), 0.0, col)
      # Express col as a linear combination of the pivot columns of rref
      linear_combination = zeros(ncols)
      linear_combination[pivots] = -col[1:length(pivots)]
      linear_combination[i] = 1.0
      push!(nullvectors, linear_combination)
    end
  end
  nullspace = reduce(hcat, nullvectors; init = Matrix{Float64}(undef, ncols, 0))

  return pivots, nullspace
end

function independent_cols!(M::AbstractMatrix, config::SparseRREFConfig)::Tuple{Vector{Int}, AbstractMatrix}
  rref, pivots = sparse_rref_with_pivots!(sparse(M), config.atol)

  ncols = size(rref, 2)
  nullvectors = Vector{Vector{Float64}}()
  for i in 1:ncols
    if i ∉ pivots
      # Column is dependent
      col = Vector(rref[:, i])
      col .= ifelse.(isapproxzero.(col), 0.0, col)
      # Express col as a linear combination of the pivot columns of rref
      linear_combination = zeros(ncols)
      linear_combination[pivots] = -col[1:length(pivots)]
      linear_combination[i] = 1.0
      push!(nullvectors, linear_combination)
    end
  end
  nullspace = reduce(hcat, nullvectors; init = Matrix{Float64}(undef, ncols, 0))

  return pivots, nullspace
end

function independent_cols!(M::AbstractMatrix, config::MixConfig)::Tuple{Vector{Int}, AbstractMatrix}
  if issparse(M)
    # For sparse matrices use sparse RREF directly.
    independent_cols!(M, config.sparse_config)
  else
    # For dense matrices, check sparsity level.
    m, n = size(M)
    nnz = count(!iszero, M)
    sparsity = 1.0 - (nnz / (m * n))

    # If matrix is highly sparse, consider sparse methods.
    if sparsity > config.sparsity
      independent_cols!(sparse(M), config.sparse_config)
    else
      independent_cols!(M, config.rref_config)
    end
  end
end

# TODO Return nullspace in `independent_cols_qr` also
"Return vector of pivots with a QR based algorithm"
function independent_cols_qr(A::Matrix; rtol::Float64 = size(A, 2) * eps(eltype(A)))::Vector{Int64}
  # Column pivoted QR decomposition.
  QR = qr(A, ColumnNorm())
  atol = rtol * abs(QR.R[1, 1])

  # Estimate rank directly from QR.
  r = something(findfirst(i -> abs(QR.R[i, i]) <= atol, 1:size(QR.R, 1)), size(A, 2) + 1) - 1

  # Return independent column indices.
  QR.p[1:r]
end

"Sparse QR-based method for selecting independent columns."
function independent_cols_qr(A::SparseMatrixCSC; rtol::Float64 = size(A, 2) * eps(eltype(A)))::Vector{Int64}
  # SuiteSparse QR factorization
  QR = qr(A)
  R = QR.R
  m, n = size(R)

  # Compute absolute tolerance based on the max diagonal element.
  atol = rtol * maximum(abs, diag(R))

  # Extract independent columns based on tolerance.
  diag_indices = min(m, n)
  independent_positions = findall(i -> abs(R[i, i]) > atol, 1:diag_indices)

  # Map back to original column indices via the column permutation.
  return QR.pcol[independent_positions]
end

"""
Construct the augmented matrix [A|b] for a linear reconciliation problem.

The augmented matrix combines the coefficient matrix A with the constant terms vector b,
forming the standard representation of a linear system Ax = b.
Variables respect the following order: measured, unmeasured, fixed.

This method is useful for debugging.
If the system is not linear an error is returned.
"""
function augmented_matrix(rec::ReconciliationProblem)
  is_linear(rec) || throw(ArgumentError("The problem is not linear."))

  (; c, A) = equations_matrix(rec)
  hcat(copy(transpose(A)), -c)
end

"Debug consistency check: verifies that the set of variables in 'rec' matches exactly the set of variables referenced across all equations. Used after equation removal to detect variables that have become orphaned."
function check_variables(rec::ReconciliationProblem)::Bool
  vars_in_eqs = Set(Iterators.flatten(rec.eq_to_var[eid] for eid in keys(rec.equations)))

  if length(vars_in_eqs) != length(rec.variables)
    return false
  end

  for v in keys(rec.variables)
    if v ∉ vars_in_eqs
      return false
    end
  end

  return true
end

"Remove redundant equations using the RREF algorithm."
function remove_redundant_equations!(
  prerec::PresolvedReconciliationProblem;
  config::PresolveConfig = PresolveConfig(),
  debug_mode::Bool = false,
)::Result{RedundancyData, SolverError}
  rec = prerec.pre
  (; c, A, linear_term_idx, bilinear_term_idx) = equations_matrix(rec)
  col_to_eid = collect(keys(rec.equations))
  if config.scale
    _, scales = scale_columns!(A)
    c = c ./ scales
  end

  # Call rref algorithm, obtaining the pivot indices and the nullspace of the matrix.
  # TODO Pass internal struct to measure time spent.
  pivots, nullspace = independent_cols!(copy(A), config)

  # Store conflicts
  conflicting_dependencies = Dictionary{EID, Vector{EID}}()
  conflicting_residuals = Dictionary{EID, Float64}()

  # Check for conflicts
  for nullvector in eachcol(nullspace)
    # Check if the nullvector is also a nullvector for the constant terms
    reference_const = dot(c, nullvector)
    if !isapproxzero(reference_const)
      # Nullvector is not nullvector for the constant terms -> conflict
      involved_cols = findall(!isapproxzero, nullvector)
      involved_equations = col_to_eid[involved_cols]
      eid = minimum(filter(eq -> !(eq in col_to_eid[pivots]), involved_equations))
      set!(conflicting_dependencies, eid, involved_equations)
      set!(conflicting_residuals, eid, reference_const)
    end
  end

  # Check if conflicting dependencies are detected
  if !isempty(conflicting_dependencies)
    return InfeasibleEquationError(
      "The left hand side of the following equations can be expressed by other equations but the right hand side don't match: $(collect(keys(conflicting_residuals)))",
      collect(keys(conflicting_residuals)),
      collect(values(conflicting_residuals)),
      Dictionary(
        ["Reason for infeasibility", "Failure point", "Conflicting equation dependencies", "Previous replacements"],
        ["Contradicting equations", "Removal of redundant equations", conflicting_dependencies, prerec.replacements],
      ),
    )
  end

  # Only keep independent equations.
  # The indices of equations which were removed corresponds to setdiff(1:length(rec.equations), pivots),
  # but for removing equations we need the corresponding EID's.
  removed_eqs_eid = [eid for (i, eid) in enumerate(keys(rec.equations)) if i ∉ pivots]
  if !isempty(removed_eqs_eid)
    # The pivots refer to equations position in the backing array,
    # so we convert to equation id before removal.
    remove_equation!(rec, removed_eqs_eid)
  end

  if debug_mode && !check_variables(rec)
    # Typically means that there exist variables which do not point to any equation,
    # this can be removed from the system.
    return InternalError("Inconsistency error in step: remove redundant equations.")
  end

  if !isnothing(prerec.rdata)
    @info "The RedundancyData already exists and will be overwritten."
  end
  prerec.rdata = RedundancyData(pivots, c, A, removed_eqs_eid, linear_term_idx, bilinear_term_idx)
end
function remove_redundant_equations!(
  rec::ReconciliationProblem;
  config::PresolveConfig = PresolveConfig(),
  debug_mode::Bool = false,
)::Result{RedundancyData, SolverError}
  prerec = PresolvedReconciliationProblem(rec)
  remove_redundant_equations!(prerec; config, debug_mode)
end

"""
In-place presolve.
"""
function presolve!(
  prerec::PresolvedReconciliationProblem;
  config::PresolveConfig = PresolveConfig(),
  # TODO Pass struct to collect intermediate timings info.
)::Result{PresolvedReconciliationProblem, SolverError}

  # Identify variables not occuring in any equation
  update_unobservable_vars!(prerec)

  @try remove_fixed_variables!(prerec)

  # Perform trivial and bilinear substitutions, which may require various passes.
  @try trivial_and_bilinear_substitutions_recursive!(prerec)

  # Remove redundant equations.
  @try remove_redundant_equations!(prerec; config)

  prerec
end
function presolve!(
  rec::ReconciliationProblem;
  config::PresolveConfig = PresolveConfig(),
  # TODO Pass struct to collect intermediate timings info.
)::Result{PresolvedReconciliationProblem, SolverError}
  prerec = PresolvedReconciliationProblem(rec)
  @try presolve!(prerec; config)
end

"""
Apply presolve scheme.
"""
function presolve(
  prerec::PresolvedReconciliationProblem;
  config::PresolveConfig = PresolveConfig(),
)::Result{PresolvedReconciliationProblem, SolverError}
  pre = deepcopy(prerec)
  @try presolve!(pre; config)
end
function presolve(
  rec::ReconciliationProblem;
  config::PresolveConfig = PresolveConfig(),
)::Result{PresolvedReconciliationProblem, SolverError}
  pre = deepcopy(rec)
  @try presolve!(pre; config)
end

function identify_determined_variables_from_equation_matrix!(
  vars_not_unique::Dictionary{Label, UnmeasuredVariable},
  cols_fixed::Set{Int64},
  problem::ReconciliationProblem,
  evaluation_values::Dictionary{String, Float64},
  eq_matrix::AbstractMatrix,
  rows_to_check::Set{Int64},
)::Bool
  found = false
  for i in rows_to_check
    row = eq_matrix[i, :]
    # Find all nonzero entries of a row as they are the relevant ones
    nonzero_entries = findall(!isapproxzero, row)
    # Determine the variables in this equation that haven't been identified as unique yet
    nonzero_nonfixed_entries = setdiff(nonzero_entries, cols_fixed)

    # Variable uniquely defined by the equation
    if length(nonzero_nonfixed_entries) == 1
      # If all but one variable in the equation is uniquely defined, the remaining one hase to be uniquely defined as well.
      varname = collect(keys(problem.unmeasured))[only(nonzero_nonfixed_entries)]
      if haskey(vars_not_unique, varname)
        delete!(vars_not_unique, varname)
        found = true
      end
      union!(cols_fixed, Set(nonzero_nonfixed_entries)) # store columns linked with uniquely defined unmeasured variables
    else
      # The equation contains several unmeasured variables.
      # Check if all terms (they are of the form `factor * variable`) are minimal or all terms are maximal
      relevant_unmeasured = collect(values(problem.unmeasured))[nonzero_nonfixed_entries]
      relevant_factors = row[nonzero_nonfixed_entries]
      all_minimal = is_minimal(relevant_unmeasured, relevant_factors, evaluation_values)
      all_maximal = is_maximal(relevant_unmeasured, relevant_factors, evaluation_values)
      if all_minimal || all_maximal
        # Variable uniquely defined because of variable bounds
        union!(cols_fixed, Set(nonzero_nonfixed_entries)) # store columns linked with uniquely defined unmeasured variables
        for j in nonzero_nonfixed_entries
          varname = collect(keys(problem.unmeasured))[j]
          if haskey(vars_not_unique, varname)
            delete!(vars_not_unique, varname)
            found = true
          end
        end
      end
    end
  end
  return found
end


"""
Find unobservable measured and unmeasured variables from the constraint equations for a ReconciliationProblem around a given point evaluation_values.
The idea is to linearize the constraints of the ReconciliationProblem around the given values (i.e., compute the jacobian) and check the linearized system for not uniquely solvable variables using the rref algorithm.
This method runs more stable if the ReconciliationProblem has been presolved before.
"""
function find_unobservable_variables_in_constraints!(
  prerec::PresolvedReconciliationProblem,
  evaluation_values::Dictionary{String, Float64} = Dictionary{String, Float64}(),
  run_presolve_first::Bool = false,
)::Dictionary{Label, UnmeasuredVariable}
  if run_presolve_first
    presolve!(prerec)
  end
  problem = prerec.pre
  (; Jx, Jy, Jz, b) = jacobian(problem, evaluation_values) # linearize problem at evaluation_values
  # Consider the jacobian for the unmeasured variables in rref in order to see
  # which unmeasured variables are uniquely determined by measured and fixed variables
  # and which depend on other unmeasured variables.
  # Using a stricter tolerance can potentially obscure linear dependencies due to numerical errors.
  # However, this seems preferable to the erroneous identification of linear dependencies (and, consequently, unmeasurable variables).
  rref, pivots = rref_with_pivots!(deepcopy(Jy), DEFAULT_ZTOL * 1e-1)

  # Label all unmeasured variables as "not unique" first
  # Remove the uniquely defined ones in the following
  vars_not_unique = deepcopy(problem.unmeasured)
  cols_fixed = Set{Int64}()

  found = true
  checkagain_row_idx_rref = Set(axes(rref, 1))
  checkagain_row_idx_original = Set(axes(Jy, 1))
  while found
    found = false
    # Check equations in rref if they determine variables due to their bounds
    found =
      identify_determined_variables_from_equation_matrix!(
        vars_not_unique,
        cols_fixed,
        problem,
        evaluation_values,
        rref,
        checkagain_row_idx_rref,
      ) || found

    # Check original equations (i.e. Jy) if they determine variables due to their bounds
    found =
      identify_determined_variables_from_equation_matrix!(
        vars_not_unique,
        cols_fixed,
        problem,
        evaluation_values,
        Jy,
        checkagain_row_idx_original,
      ) || found

    if found
      # Identifying a variable as not free can cause other variables to be not free either.
      # Hence, all equations containing variables that are uniquely defined because of variable bounds have to be checked again.
      # Indices of rows in which a variable of cols_fixed participates:
      checkagain_row_idx_rref = Set{Int64}()
      checkagain_row_idx_original = Set{Int64}()
      for colidx in cols_fixed
        union!(checkagain_row_idx_rref, Set(findall(!isapproxzero, rref[:, colidx])))
        union!(checkagain_row_idx_original, Set(findall(!isapproxzero, Jy[:, colidx])))
      end
    end
  end

  merge!(prerec.unobservable, vars_not_unique)
end
function find_unobservable_variables_in_constraints(
  prerec::PresolvedReconciliationProblem,
  evaluation_values::Dictionary{String, Float64} = Dictionary{String, Float64}(),
  run_presolve_first::Bool = false,
)::Dictionary{Label, UnmeasuredVariable}
  pre = deepcopy(prerec)
  find_unobservable_variables_in_constraints!(pre, evaluation_values, run_presolve_first)
end
function find_unobservable_variables_in_constraints(
  rec::ReconciliationProblem,
  evaluation_values::Dictionary{String, Float64} = Dictionary{String, Float64}(),
  run_presolve_first::Bool = false,
)::Dictionary{Label, UnmeasuredVariable}
  prerec = PresolvedReconciliationProblem(rec)
  find_unobservable_variables_in_constraints(prerec, evaluation_values, run_presolve_first)
end

"""
Find free measured and unmeasured variables from the constraint equations for a ReconciliationSolution around a computed solution values.
"""
function find_unobservable_variables_in_constraints!(rsol::ReconciliationSolution; run_presolve_first::Bool = false)
  prerec = PresolvedReconciliationProblem(rsol.problem; unobservable = rsol.unobservable_variables)
  find_unobservable_variables_in_constraints!(prerec, rsol.variables_value, run_presolve_first)
end
function find_unobservable_variables_in_constraints(rsol::ReconciliationSolution; run_presolve_first::Bool = false)
  prerec = PresolvedReconciliationProblem(rsol.problem; unobservable = rsol.unobservable_variables)
  find_unobservable_variables_in_constraints(prerec, rsol.variables_value, run_presolve_first)
end

end
