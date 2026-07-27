"""
Variables and equations which are used to represent model constraints, typically:

- Mass-balance equations.
- Transfer coefficients.
- Coupling between different substance layers.
- Coupling between different time periods.
"""
module ConstraintsModule

using Dictionaries
using StructEquality
using ..ComparisonModule
using ..ErrorModule

export
  Variable,
  Label,
  UnmeasuredVariable,
  MeasuredVariable,
  FixedVariable,
  Term,
  ConstantTerm,
  LinearTerm,
  BilinearTerm,
  Equation,
  is_linear,
  is_bilinear,
  lower_bound,
  upper_bound,
  value,
  get_index,
  is_transfer_coefficient,
  satisfies_upper_bound,
  satisfies_lower_bound,
  satisfies_bounds,
  matches_upper_bound,
  matches_lower_bound,
  matches_bound,
  is_minimal,
  is_maximal,
  assign_mean!,
  delete_unobservable_variable_values!,
  term_names,
  term_names!,
  evaluate

# ==============================
# Different types of Variables
# ==============================

"Label used for variable names."
const Label = String

"""
Variables with their respective uncertainties, if available.
"""
abstract type Variable end

"""
Variable that is measured with an associated uncertainty.
Normality assumption is impicit.
"""
@struct_hash_equal struct MeasuredVariable <: Variable
  "Name of this variable."
  name::Label
  "Value associated to the variable (typically mean value of a normal distribution)."
  value::Float64
  "Uncertainty associated to the variable (typically standard deviation of a normal distribution)."
  uncertainty::Float64
  "Additional information on this variable (typically name in extended form)."
  info::String
  "Physical units (if available)."
  units::Union{String, Missing}
  "Whether the variable represents a transfer coefficient (if available)."
  tc::Union{Bool, Missing}
  "Admissible lower bound (if available)."
  lb::Union{Float64, Missing}
  "Admissible upper bound (if available)."
  ub::Union{Float64, Missing}
end

function MeasuredVariable(;
  name::Label,
  value::Float64,
  uncertainty::Float64,
  info::String = "",
  units::Union{String, Missing} = missing,
  tc::Union{Bool, Missing} = missing,
  lb::Union{Float64, Missing} = missing,
  ub::Union{Float64, Missing} = missing)
  # By default:
  # - Transfer coefficients are only allowed in the [0, 1] interval.
  # - In general, if no lower bound is specified, we assume it is a flow with nonnegative value.
  tc = ismissing(tc) ? occursin("Transfer coefficient", info) : tc
  if tc
    lb = ismissing(lb) ? zero(Float64) : lb
    ub = ismissing(ub) ? one(Float64) : ub
  else
    if ismissing(lb)
      lb = zero(Float64)
    end
  end
  if !ismissing(lb)
    if !relaxed_leq(lb, value)
      throw(
        IllposedProblemError(
          "Measured value cannot be smaller than the lower bound for variable $name.",
          "",
          Dictionary{String, Any}(),
        ),
      )
    elseif !ismissing(ub)
      if !relaxed_geq(ub, value)
        throw(
          IllposedProblemError(
            "Measured value cannot be greater than the upper bound for variable $name.",
            "",
            Dictionary{String, Any}(),
          ),
        )
      elseif relaxed_isapprox(lb, ub)
        @debug "Made variable $name fixed as lower and upper bound conincide."
        return FixedVariable(; name, value, info, units)
      end
    end
  end
  MeasuredVariable(name, value, uncertainty, info, units, tc, lb, ub)
end

value(var::MeasuredVariable) = var.value
is_transfer_coefficient(var::MeasuredVariable) = var.tc
lower_bound(var::MeasuredVariable)::Union{Float64, Missing} = var.lb
upper_bound(var::MeasuredVariable)::Union{Float64, Missing} = var.ub

"Assign measured variables to their mean."
function assign_mean!(
  measured::Dictionary{Label, MeasuredVariable},
  values::Dictionary{Label, Float64} = Dictionary{Label, Float64}(),
)::Dictionary{Label, Float64}
  for var_name in keys(measured)
    if !haskey(values, var_name)
      var = measured[var_name]
      mean = var.value
      set!(values, var_name, mean)
    end
  end
  values
end

"""
Variable for which there is no data available.
"""
@struct_hash_equal struct UnmeasuredVariable <: Variable
  "Name of this variable."
  name::Label
  "Additional information on this variable (typically name in extended form)."
  info::String
  "Physical units (if available)."
  units::Union{String, Missing}
  "Whether the variable represents a transfer coefficient (if available)."
  tc::Union{Bool, Missing}
  "Admissible lower bound (if available)."
  lb::Union{Float64, Missing}
  "Admissible upper bound (if available)."
  ub::Union{Float64, Missing}
end

function UnmeasuredVariable(;
  name::Label,
  info::String = "",
  units::Union{String, Missing} = missing,
  tc::Union{Bool, Missing} = missing,
  lb::Union{Float64, Missing} = missing,
  ub::Union{Float64, Missing} = missing,
)
  # By default:
  # - Transfer coefficients are only allowed in the [0, 1] interval.
  # - In general, if no lower bound is specified, we assume it is a flow with nonnegative value.
  tc = ismissing(tc) ? occursin("Transfer coefficient", info) : tc
  if tc
    lb = ismissing(lb) ? zero(Float64) : lb
    ub = ismissing(ub) ? one(Float64) : ub
  else
    if ismissing(lb)
      lb = zero(Float64)
    end
  end
  if !ismissing(lb) && !ismissing(ub)
    if !relaxed_leq(lb, ub)
      throw(
        IllposedProblemError(
          "Lower bound cannot be greater than upper bound for variable $name.",
          "",
          Dictionary{String, Any}(),
        ),
      )
    elseif relaxed_isapprox(lb, ub)
      @debug "Made variable $name fixed as lower and upper bound conincide."
      return FixedVariable(; name, value = lb, info, units)
    end
  end
  UnmeasuredVariable(name, info, units, tc, lb, ub)
end

value(::UnmeasuredVariable) = nothing
is_transfer_coefficient(var::UnmeasuredVariable) = var.tc
lower_bound(var::UnmeasuredVariable)::Union{Float64, Missing} = var.lb
upper_bound(var::UnmeasuredVariable)::Union{Float64, Missing} = var.ub

"Remove values of unobservable variables from the solution dictionary."
function delete_unobservable_variable_values!(
  variables_value::Dictionary{Label, Float64},
  unobservable_variables::Dictionary{Label, UnmeasuredVariable},
)
  for var in collect(keys(unobservable_variables))
    if haskey(variables_value, var)
      delete!(variables_value, var)
    end
  end
end

"""
Variable (or parameter) for which there is a certain and fixed value.
"""
@struct_hash_equal @kwdef struct FixedVariable <: Variable
  "Name of this variable."
  name::Label
  "Value assigned to the variable."
  value::Float64
  "Additional information on this variable (typically name in extended form)."
  info::String = ""
  "Physical units (if available)."
  units::Union{String, Missing} = ""
  "Whether the variable represents a transfer coefficient (if available)."
  tc::Bool = false
end
FixedVariable(name::Label, value::Float64) = FixedVariable(; name, value, info = "", units = missing)

value(var::FixedVariable) = var.value
is_transfer_coefficient(var::FixedVariable) = var.tc
upper_bound(var::FixedVariable) = var.value
lower_bound(var::FixedVariable) = var.value

"Check if a value satisfies the upper bound of a variable (if given)."
function satisfies_upper_bound(var::Variable, value::Float64)::Bool
  if ismissing(upper_bound(var)) || relaxed_leq(value, upper_bound(var))
    return true
  else
    @debug "Upper bound $(upper_bound(var)) for $(typeof(var)) variable $(var.name) not satisfied, value $(value)."
    return false
  end
end

"Check if a value satisfies the lower bound of a variable (if given)."
function satisfies_lower_bound(var::Variable, value::Float64)::Bool
  if ismissing(lower_bound(var)) || relaxed_geq(value, lower_bound(var))
    return true
  else
    @debug "Lower bound $(lower_bound(var)) for $(typeof(var)) variable $(var.name) not satisfied, value $(value)."
    return false
  end
end

"Check if a value satisfies the upper and lower bounds of a variable (if given)."
function satisfies_bounds(var::Variable, value::Float64)::Bool
  return satisfies_upper_bound(var, value) && satisfies_lower_bound(var, value)
end

"Check if a value matches the upper bound of a variable (if given)."
function matches_upper_bound(var::Variable, value::Float64)::Bool
  if ismissing(upper_bound(var))
    return false
  else
    return relaxed_isapprox(value, upper_bound(var))
  end
end

function matches_upper_bound(vars::Vector{<:Variable}, values::Dictionary{Label, Float64})::Bool
  # Iterate over all variables
  for var in vars
    # If a variable doesn't match the upper bound, return false
    if !haskey(values, var.name) || !matches_upper_bound(var, values[var.name])
      return false
    end
  end
  # As no variable was found that doesn't match the upper bound, return true
  return true
end

"Check if a value matches the lower bound of a variable (if given)."
function matches_lower_bound(var::Variable, value::Float64)::Bool
  if ismissing(lower_bound(var))
    return false
  else
    return relaxed_isapprox(value, lower_bound(var))
  end
end

function matches_lower_bound(vars::Vector{<:Variable}, values::Dictionary{Label, Float64})::Bool
  # Iterate over all variables
  for var in vars
    # If a variable doesn't match the lower bound, return false
    if !haskey(values, var.name) || !matches_lower_bound(var, values[var.name])
      return false
    end
  end
  # As no variable was found that doesn't match the lower bound, return true
  return true
end

"Check if a value matches the upper or the lower bound of a variable (if given)."
function matches_bound(var::Variable, value::Float64)::Bool
  return matches_upper_bound(var, value) || matches_lower_bound(var, value)
end

function matches_bound(vars::Vector{<:Variable}, values::Dictionary{Label, Float64})::Bool
  # Iterate over all variables
  for var in vars
    # If a variable doesn't match a bound, return false
    if !haskey(values, var.name) || !matches_bound(var, values[var.name])
      return false
    end
  end
  # As no variable was found that doesn't match a bound, return true
  return true
end

"Check if `factor * value` is the minimal value that `factor * var` can attain (due to variable bounds)."
function is_minimal(var::Variable, factor::Float64, value::Float64)::Bool
  if iszero(factor)
    return true
  elseif factor > 0
    return matches_lower_bound(var, value)
  else
    # factor < 0
    return matches_upper_bound(var, value)
  end
end

function is_minimal(vars::Vector{<:Variable}, factors::Vector{Float64}, values::Dictionary{Label, Float64})::Bool
  nvars = length(vars)
  @assert length(factors) == nvars

  # Iterate over all variables
  for i in 1:nvars
    var = vars[i]
    factor = factors[i]
    # If a variable is not minimal, return false
    if !haskey(values, var.name) || !is_minimal(var, factor, values[var.name])
      return false
    end
  end
  # As no variable was found that isn't minimal, return true
  return true
end

"Check if `factor * value` is the maximal value that `factor * var` can attain (due to variable bounds)."
function is_maximal(var::Variable, factor::Float64, value::Float64)::Bool
  if iszero(factor)
    return true
  elseif factor > 0
    return matches_upper_bound(var, value)
  else
    # factor < 0
    return matches_lower_bound(var, value)
  end
end

function is_maximal(vars::Vector{<:Variable}, factors::Vector{Float64}, values::Dictionary{Label, Float64})::Bool
  nvars = length(vars)
  @assert length(factors) == nvars

  # Iterate over all variables
  for i in 1:nvars
    var = vars[i]
    factor = factors[i]
    # If a variable is not maximal, return false
    if !haskey(values, var.name) || !is_maximal(var, factor, values[var.name])
      return false
    end
  end
  # As no variable was found that isn't maximal, return true
  return true
end

# ====================
# Terms and equations
# ====================

"""
Represents a term (typically a multivariate monomial) as part of an equation. 
"""
abstract type Term end

"A constant term is one of the form `factor`. The name is optional."
@struct_hash_equal struct ConstantTerm <: Term
  "Constant value."
  value::Float64
end
ConstantTerm() = ConstantTerm(zero(Float64))

function Base.:(+)(c1::ConstantTerm, c2::ConstantTerm)::ConstantTerm
  ConstantTerm(c1.value + c2.value)
end
function Base.iszero(c1::ConstantTerm)
  iszero(c1.value)
end

"A linear term is one of the form `factor * variable`."
@struct_hash_equal struct LinearTerm <: Term
  "Name of the variable."
  name::Label
  "Multiplicative factor (coefficient)."
  factor::Float64
  function LinearTerm(name, factor::Float64)
    if iszero(factor)
      throw(ArgumentError("Factor cannot be zero."))
    end
    new(name, factor)
  end
end

"A bilinear term is one of the form `factor * variable1 * variable2`."
@struct_hash_equal struct BilinearTerm <: Term
  "Name of the first variable."
  name1::Label
  "Name of the second variable."
  name2::Label
  "Multiplicative factor (coefficient)."
  factor::Float64
  function BilinearTerm(name1, name2, factor::Float64)
    if iszero(factor)
      throw(ArgumentError("Factor cannot be zero."))
    end
    new(name1, name2, factor)
  end
end

"""
Return index of a LinearTerm with variable 'var' within a list of LinearTerms 'terms'.
Terminate when the first LinearTerm with variable 'var' is found
(i.e. assume that there aren't several linear terms with the same variable).
Return 'nothing' if there is no LinearTerm with variable 'var'.
"""
function get_index(terms::Vector{LinearTerm}, var::Label)::Union{Int64, Nothing}
  idx = nothing
  for (i, lt) in enumerate(terms)
    if lt.name == var
      idx = i
      break
    end
  end
  return idx
end

"""
Return index of a BilinearTerm with product variables 'var1' and 'var2' within a list of BilinearTerms 'terms'.
Terminate when the first BilinearTerm with variables 'var1' and 'var2' (i.e. 'var1 * var2' or 'var2 * var1') is found
(i.e. assume that there aren't several bilinear terms with the same variables).
Return 'nothing' if there is no BilinearTerm with variables 'var1' and 'var2'.
"""
function get_index(terms::Vector{BilinearTerm}, var1::Label, var2::Label)::Union{Int64, Nothing}
  idx = nothing
  for (i, bt) in enumerate(terms)
    if (bt.name1 == var1 && bt.name2 == var2) || (bt.name1 == var2 && bt.name2 == var1)
      idx = i
      break
    end
  end
  return idx
end

"""
A (possibly nonlinear) equation of the form `C + ∑ᵢ L[i] + ∑ᵢ B[i] = 0` where:
- `C[i]` are constant terms
- `L[i]` are linear terms
- `B[i]` are nonlinear terms
"""
@struct_hash_equal mutable struct Equation
  "Constant term of the equation."
  constant_term::ConstantTerm
  "Linear terms of the equation."
  const linear_terms::Vector{LinearTerm}
  "Bilinear (nonlinear) terms of the equation."
  const bilinear_terms::Vector{BilinearTerm}
end

function Equation(linear_terms::Vector{LinearTerm}, bilinear_terms::Vector{BilinearTerm})
  Equation(ConstantTerm(), linear_terms, bilinear_terms)
end

function is_linear(eq::Equation)::Bool
  isempty(eq.bilinear_terms)
end
function is_bilinear(eq::Equation)::Bool
  !is_linear(eq)
end

"Return the label for each named term in the equation."
function term_names(eq::Equation)::Vector{Label}
  labels = Label[]
  term_names!(labels, eq)
end
function term_names!(labels::Vector{Label}, eq::Equation)::Vector{Label}
  for t in eq.linear_terms
    push!(labels, t.name)
  end
  for t in eq.bilinear_terms
    push!(labels, t.name1)
    push!(labels, t.name2)
  end
  unique!(labels)
end
"Return the label for each named term in the vector of linear terms."
function term_names(linear_terms::Vector{LinearTerm})::Vector{Label}
  term_names(Equation(linear_terms, BilinearTerm[]))
end

function Base.show(io::IO, eq::Equation)
  parts = String[]

  # Add linear terms.
  for term in eq.linear_terms
    coeff = term.factor
    var = term.name

    p = if coeff == 1.0
      " + $(var)"
    elseif coeff == -1.0
      " - $(var)"
    elseif coeff > 0
      " + $(coeff) * $(var)"
    else
      "$(coeff) * $(var)"
    end
    push!(parts, p)
  end

  # Add bilinear terms.
  for term in eq.bilinear_terms
    coeff = term.factor
    p = if coeff == 1.0
      " + $(term.name1) * $(term.name2)"
    elseif coeff == -1.0
      " - $(term.name1) * $(term.name2)"
    elseif coeff > 0
      " + $(coeff) * $(term.name1) * $(term.name2)"
    else
      "$(coeff) * $(term.name1) * $(term.name2)"
    end
    if isempty(parts) && coeff >= 0
      p *= " + "
    end
    push!(parts, p)
  end

  # Equations without linear and bilinear terms (i.e., lhs is '0')
  if isempty(parts)
    push!(parts, "0")
  end

  # Add constant term.
  val = eq.constant_term.value
  p = if val > 0
    " = -$val"
  elseif iszero(val)
    " = 0"
  else
    " = $(abs(val))"
  end
  push!(parts, p)

  # The left-most term doesn't need a `+` symbol.
  first_term = first(parts)
  if first(strip(first_term)) == '+'
    parts[1] = split(first_term, '+')[2]
  end
  equation_str = join(parts)
  print(io, equation_str)
end

# ======================
# Numerical evaluation
# ======================

"Evaluate the constant term."
function evaluate(t::ConstantTerm)::Float64
  t.value
end

"Evaluate the linear term."
function evaluate(t::LinearTerm, var::Variable)::Float64
  if t.name != var.name
    @warn "The name of the term $(t.name) and the name of the variable $(var.name) don't match."
  end
  t.factor * var.value
end
function evaluate(t::LinearTerm, vars::Vector{<:Variable})::Float64
  idx = findall(v -> v.name == t.name, vars)
  if isnothing(idx)
    throw(ArgumentError("Variable $(t.name) not found in the result."))
  end
  evaluate(t, vars[only(idx)])
end
function evaluate(t::LinearTerm, vars::Dictionary{Label, Float64})::Float64
  t.factor * vars[t.name]
end

"Evaluate the bilinear term."
function evaluate(t::BilinearTerm, var1::Variable, var2::Variable)::Float64
  if t.name1 != var1.name
    @warn "The name of the term $(t.name1) and the name of the variable $(var1.name) don't match."
  end
  if t.name2 != var2.name
    @warn "The name of the term $(t.name2) and the name of the variable $(var2.name) don't match."
  end
  t.factor * var1.value * var2.value
end
function evaluate(t::BilinearTerm, vars::Vector{<:Variable})::Float64
  idx1 = findall(v -> v.name == t.name1, vars)
  if isnothing(idx1)
    throw(ArgumentError("Variable $(t.name1) not found in the result."))
  end
  idx2 = findall(v -> v.name == t.name2, vars)
  if isnothing(idx2)
    throw(ArgumentError("Variable $(t.name2) not found in the result."))
  end
  evaluate(t, vars[only(idx1)], vars[only(idx2)])
end
function evaluate(t::BilinearTerm, vars::Dictionary{Label, Float64})::Float64
  val1 = get(vars, t.name1, nothing)
  val2 = get(vars, t.name2, nothing)
  # Cover for cases in which some variable is missing from the assignment.
  # Typically when unobservable variables are present.
  if !isnothing(val1) && !isnothing(val2)
    t.factor * val1 * val2
  elseif isnothing(val1) && !isnothing(val2) && iszero(val2)
    0.0
  elseif !isnothing(val1) && iszero(val1) && isnothing(val2)
    0.0
  else
    throw(ArgumentError("Result cannot be determined since $(t.name1) and $(t.name2) not found in the result."))
  end
end

"Evaluate a term consisting of a ConstantTerm and a vector of LinearTerms."
function evaluate(constant_term::ConstantTerm, linear_terms::Vector{LinearTerm}, vars::Vector{<:Variable})::Float64
  acc = zero(Float64)

  # Accumulate over constant terms.
  acc += evaluate(constant_term)

  # Accumulate over linear terms.
  for t in linear_terms
    acc += evaluate(t, vars)
  end
  acc
end
function evaluate(
  constant_term::ConstantTerm,
  linear_terms::Vector{LinearTerm},
  vars::Dictionary{Label, Float64},
)::Float64
  acc = zero(Float64)

  # Accumulate over constant terms.
  acc += evaluate(constant_term)

  # Accumulate over linear terms.
  for t in linear_terms
    acc += evaluate(t, vars)
  end
  acc
end

"Evaluate the equation."
function evaluate(eq::Equation, vars::Vector{<:Variable})::Float64
  acc = zero(Float64)

  # Accumulate over constant terms.
  acc += evaluate(eq.constant_term)

  # Accumulate over linear terms.
  for t in eq.linear_terms
    acc += evaluate(t, vars)
  end

  # Accumulate over bilinear terms.
  for t in eq.bilinear_terms
    acc += evaluate(t, vars)
  end
  acc
end
function evaluate(eq::Equation, vars::Dictionary{Label, Float64})::Float64
  acc = zero(Float64)

  # Accumulate over constant terms.
  acc += evaluate(eq.constant_term)

  # Accumulate over linear terms.
  for t in eq.linear_terms
    acc += evaluate(t, vars)
  end

  # Accumulate over bilinear terms.
  for t in eq.bilinear_terms
    acc += evaluate(t, vars)
  end
  acc
end

end
