"""
Parse equations and bound constraints from string format into a more structured format.
"""
module ConstraintParserModule

using ..ComparisonModule
using ..ConstraintsModule

export
  determine_constraint_type,
  parse_constraint,
  parse_equation,
  parse_bound_constraint,
  BoundConstraint

"Simplified parser with unified equation handling."
abstract type ParsedConstraint end

"Type for bound constraints (lower <= variable <= upper)."
struct BoundConstraint <: ParsedConstraint
  name::String
  lower_bound::Union{Float64, Missing}
  upper_bound::Union{Float64, Missing}
end

"Helper function to determine the constraint type."
function determine_constraint_type(constraint::String)::Union{Type{BoundConstraint}, Type{Equation}}
  # Remove all spaces for easier parsing.
  clean_constraint = replace(constraint, r"\s+" => "")

  if occursin("<=", clean_constraint) || occursin("≤", clean_constraint) || occursin("<", clean_constraint) ||
     occursin(">=", clean_constraint) || occursin("≥", clean_constraint) || occursin(">", clean_constraint)
    BoundConstraint
  elseif occursin("=", clean_constraint) && !occursin("<=", clean_constraint) && !occursin("≤", clean_constraint) &&
         !occursin(">=", clean_constraint) && !occursin("≥", clean_constraint)
    Equation
  else
    throw(ArgumentError("Unrecognized constraint format: $constraint"))
  end
end

"Main parsing function that handles all constraint types."
function parse_constraint(constraint::String)::Union{Equation, BoundConstraint}
  constraint_type = determine_constraint_type(constraint)

  if constraint_type == BoundConstraint
    parse_bound_constraint(constraint)
  elseif constraint_type == Equation
    parse_equation(constraint)
  else
    throw(ArgumentError("Unknown constraint type"))
  end
end

"Helper to extract the coefficient and variables."
function parse_term(term::String)::Term
  factor = one(Float64)
  sign = one(Float64)
  # Check if the term starts with `+` or `-`, adjust the sign accordingly
  # and then remove it from the term if present.
  if startswith(term, '-')
    sign = -one(Float64)
    term = term[2:end]
  elseif startswith(term, '+')
    term = term[2:end]
  end

  # Split according to `*` into components.
  components = split(term, '*')

  factor = sign * one(Float64)
  vars_components = []

  for comp in components
    c = tryparse(Float64, comp)
    if !isnothing(c)
      factor = factor * c
    else
      # Split according to `/` into components.
      division_components = split(comp, '/')
      # Parse numerator
      c = tryparse(Float64, division_components[1])
      if !isnothing(c)
        factor = factor * c
      else
        push!(vars_components, division_components[1])
      end
      # Parse denominator(s)
      for denom in division_components[2:end]
        d = tryparse(Float64, denom)
        if isnothing(d)
          # Allow numbers only as denominator
          throw(ArgumentError("Unsupported term structure: $term"))
        end
        factor = factor / d
      end
    end
  end

  result = if length(vars_components) == 0
    # Constant term (e.g., 2).
    # Note: These terms are not typically present in STAN models, 
    # but are kept here to allow broader use of this method.
    ConstantTerm(factor)
  elseif length(vars_components) == 1
    # Linear term (e.g., FV:726 or 3 * FV:726).
    if isapproxzero(factor)
      ConstantTerm(0.0)
    else
      LinearTerm(vars_components[1], factor)
    end
  elseif length(vars_components) == 2
    # Bilinear term
    if isapproxzero(factor)
      ConstantTerm(0.0)
    else
      BilinearTerm(vars_components[1], vars_components[2], factor)
    end
  else
    throw(ArgumentError("Unrecognized term structure: $term"))
  end

  result
end

"Helper function to parse expression into terms."
function parse_expression(expr::AbstractString)
  # Split the expression by `+` and `-`, keeping signs attached to terms.
  terms = eachmatch(r"[+-]?[^+-]+", expr)
  terms = String.(getfield.(collect(terms), :match))

  # Initialize vectors for constant, linear and bilinear terms.
  constant_term = ConstantTerm()
  linear_terms = LinearTerm[]
  bilinear_terms = BilinearTerm[]

  # Process each term and categorize it.
  for term in terms
    parsed_term = parse_term(term)
    if isa(parsed_term, ConstantTerm)
      constant_term += parsed_term
    elseif isa(parsed_term, LinearTerm)
      push!(linear_terms, parsed_term)
    elseif isa(parsed_term, BilinearTerm)
      push!(bilinear_terms, parsed_term)
    else
      throw(ArgumentError("Unrecognized term."))
    end
  end

  return constant_term, linear_terms, bilinear_terms
end

"Parse any equality constraint (assignments, equations) into unified form: LHS - RHS = 0."
function parse_equation(constraint::String)::Equation
  # Remove spaces and split into left and right parts.
  constraint = replace(constraint, r"\s+" => "")
  left, right = split(constraint, "=")

  # Parse both sides.
  left_const, left_linear, left_bilinear = parse_expression(left)
  right_const, right_linear, right_bilinear = parse_expression(right)

  # Move everything to left side: LHS - RHS = 0
  # Constant terms: left_const - right_const
  final_constant = ConstantTerm(left_const.value - right_const.value)

  # Linear terms: combine and negate right side terms
  final_linear = copy(left_linear)
  for term in right_linear
    push!(final_linear, LinearTerm(term.name, -term.factor))
  end

  # Bilinear terms: combine and negate right side terms  
  final_bilinear = copy(left_bilinear)
  for term in right_bilinear
    push!(final_bilinear, BilinearTerm(term.name1, term.name2, -term.factor))
  end

  Equation(final_constant, final_linear, final_bilinear)
end

"Parse bound constraints with single or double inequality symbols."
function parse_bound_constraint(constraint::AbstractString)::BoundConstraint
  # Remove spaces
  constraint = replace(constraint, r"\s+" => "")

  # Handle double inequalities: lower op1 variable op2 upper
  # Count inequality operators
  inequality_count = count(c -> c in "<≤≥>", constraint)

  if inequality_count == 2
    # Pattern for double inequality: number op variable op number
    # Since spaces are removed, no \s* needed
    pattern = r"([+-]?[\d.]+)(<=|>=|≤|≥|<|>)([a-zA-Z_][a-zA-Z0-9_:]*)(<=|>=|≤|≥|<|>)([+-]?[\d.]+)"
    m = match(pattern, constraint)

    if !isnothing(m)
      lower_val = parse(Float64, m.captures[1])
      lower_op = m.captures[2]
      variable = m.captures[3]
      upper_op = m.captures[4]
      upper_val = parse(Float64, m.captures[5])

      # Validate that this is a proper double bound
      # lower_op should be <= or < (pointing right)
      # upper_op should be <= or < (pointing right)
      if lower_op in ["<", "<=", "≤"] && upper_op in ["<", "<=", "≤"]
        # Standard form: lower <= variable <= upper
        return BoundConstraint(variable, lower_val, upper_val)
      else
        throw(ArgumentError("Invalid double inequality format: $constraint"))
      end
    else
      throw(ArgumentError("Could not parse double inequality: $constraint"))
    end
  elseif inequality_count == 1
    # Handle single inequalities
    return parse_single_inequality(constraint)
  else
    throw(ArgumentError("Invalid number of inequality operators in constraint: $constraint"))
  end
end

function parse_single_inequality(constraint::AbstractString)::BoundConstraint
  # Determine which operator we have and split accordingly
  if occursin(">=", constraint)
    parts = split(constraint, ">=", limit = 2)
    operator = ">="
  elseif occursin("≥", constraint)
    parts = split(constraint, "≥", limit = 2)
    operator = "≥"
  elseif occursin("<=", constraint)
    parts = split(constraint, "<=", limit = 2)
    operator = "<="
  elseif occursin("≤", constraint)
    parts = split(constraint, "≤", limit = 2)
    operator = "≤"
  elseif occursin(">", constraint)
    parts = split(constraint, ">", limit = 2)
    operator = ">"
  elseif occursin("<", constraint)
    parts = split(constraint, "<", limit = 2)
    operator = "<"
  else
    throw(ArgumentError("No inequality operator found in constraint: $constraint"))
  end

  if length(parts) != 2
    throw(ArgumentError("Invalid inequality constraint format: $constraint"))
  end

  left, right = parts
  left_is_num = !isnothing(tryparse(Float64, left))
  right_is_num = !isnothing(tryparse(Float64, right))

  # Determine variable and bound based on which side is numeric
  if left_is_num && !right_is_num
    # Pattern: number op variable
    bound_val = parse(Float64, left)
    variable = right

    if operator in ["<=", "<", "≤"]
      # bound <= variable  →  variable >= bound (lower bound)
      BoundConstraint(variable, bound_val, missing)
    else  # operator in [">=", ">"]
      # bound >= variable  →  variable <= bound (upper bound)
      BoundConstraint(variable, missing, bound_val)
    end

  elseif !left_is_num && right_is_num
    # Pattern: variable op number
    variable = left
    bound_val = parse(Float64, right)

    if operator in ["<=", "<", "≤"]
      # variable <= bound (upper bound)
      BoundConstraint(variable, missing, bound_val)
    else  # operator in [">=", ">"]
      # variable >= bound (lower bound)
      BoundConstraint(variable, bound_val, missing)
    end

  else
    throw(
      ArgumentError("Invalid constraint format - both sides are $(left_is_num ? "numeric" : "variables"): $constraint"),
    )
  end
end

end
