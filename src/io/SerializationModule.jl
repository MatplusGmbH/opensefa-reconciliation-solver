"""
Methods to load reconciliation problems from JSON files and export solutions.

- Parse models described in JSON to a ReconciliationProblem.
- Serialize solutions to JSON.
"""
module SerializationModule

import JSON3
import Pkg
using Accessors: @set
using Dictionaries

using ..ComparisonModule
using ..ConstraintsModule
using ..ErrorInfosModule
using ..ConstraintParserModule
using ..DefaultSolverModule
using ..JuMPInterfaceModule
using ..ReconciliationModule

export load_model, serialize, load_solver, to_dict

function parse_bounds(jvar)::Tuple{Union{Float64, Missing}, Union{Float64, Missing}}
  # Parse bounds if they exist in the JSON, otherwise return missing to let the constructor apply defaults
  lb = if haskey(jvar, "lb") && !isnothing(jvar["lb"])
    jvar["lb"] == "-Inf" ? -Inf : convert(Float64, jvar["lb"])
  else
    missing
  end

  ub = if haskey(jvar, "ub") && !isnothing(jvar["ub"])
    jvar["ub"] in ("Inf", "+Inf") ? Inf : convert(Float64, jvar["ub"])
  else
    missing
  end

  lb, ub
end

function load_model(json_data::JSON3.Object)::AbstractReconciliationProblem
  name = json_data["name"]
  variables = json_data["variables"]
  equations = json_data["equations"]

  rec = ReconciliationProblem(; name)

  for jvar in variables
    type = jvar["type"]
    name = jvar["name"]
    # If no bounds are defined, use the same defaults as the variable constructors.
    var = if type == "MeasuredVariable"
      lb, ub = parse_bounds(jvar)
      MeasuredVariable(;
        name,
        value = convert(Float64, jvar["value"]),
        uncertainty = convert(Float64, jvar["uncertainty"]),
        tc = get(jvar, "is_transfer_coefficient", missing),
        lb,
        ub,
      )
    elseif type == "UnmeasuredVariable"
      lb, ub = parse_bounds(jvar)
      UnmeasuredVariable(;
        name,
        tc = get(jvar, "is_transfer_coefficient", missing),
        lb,
        ub,
      )
    elseif type == "FixedVariable"
      FixedVariable(;
        name,
        value = convert(Float64, jvar["value"]),
        tc = jvar["is_transfer_coefficient"])
    else
      throw(ArgumentError("Unexpected type $type."))
    end
    add_variable!(rec, var)
  end

  for jeq in equations
    constant_term = ConstantTerm(convert(Float64, jeq["constant_term"]["value"]))
    linear_terms = [
      LinearTerm(t["name"], convert(Float64, t["factor"])) for
      t in jeq["linear_terms"] if !isapproxzero(convert(Float64, t["factor"]))
    ]
    bilinear_terms =
      [
        BilinearTerm(t["name1"], t["name2"], convert(Float64, t["factor"])) for
        t in jeq["bilinear_terms"] if !isapproxzero(convert(Float64, t["factor"]))
      ]
    eq = Equation(constant_term, linear_terms, bilinear_terms)
    add_equation!(rec, eq; strict = true)
  end

  if haskey(json_data, "constraints")
    for constraint in json_data["constraints"]
      found = parse_constraint(constraint)
      if found isa Equation
        # For each term name, add it to the problem if it is not present.
        for name in term_names(found)
          if name ∉ keys(rec.variables)
            add_variable!(rec, UnmeasuredVariable(; name))
          end
        end
        add_equation!(rec, found)
      elseif found isa BoundConstraint
        # If the constraint defines a bound, add a new variable with such bounds to the system.
        # Also handle the case in which the variable already exists in the system and this constraint defines additional bounds.
        name = found.name
        if name ∉ keys(rec.variables)
          add_variable!(rec, UnmeasuredVariable(; name, lb = found.lower_bound, ub = found.upper_bound))
        else
          # Update bounds if this variable is already present in the model.
          var = get_variable(rec, name)
          if var isa MeasuredVariable
            var = @set var.lb = found.lower_bound
            var = @set var.ub = found.upper_bound
            rec.measured[name] = var
          elseif var isa UnmeasuredVariable
            var = @set var.lb = found.lower_bound
            var = @set var.ub = found.upper_bound
            rec.unmeasured[name] = var
          else
            throw(ArgumentError("Type of variable not supported to change bounds"))
          end
        end
      else
        throw(ArgumentError("Type of constraint not supported"))
      end
    end
  end

  rec
end
function load_model(json_data::Dict)::AbstractReconciliationProblem
  load_model(JSON3.read(JSON3.write(json_data)))
end

"Load model from JSON file."
function load_model(file::String)::AbstractReconciliationProblem
  if !isfile(file)
    throw(ArgumentError("File $file not found."))
  end
  json_data = JSON3.read(file)

  # TODO Validate all fields.
  if haskey(json_data, :content)
    # Models generated from the editor need to extract from the content field.
    json_data = json_data[:content]
  end
  load_model(json_data)
end

"If the JSON model specifies a solver, use it. Otherwise pass a default solver."
function load_solver(file::String)::ReconciliationSolver
  if !isfile(file)
    throw(ArgumentError("File $file not found."))
  end
  json_data = JSON3.read(file)

  if haskey(json_data, :content)
    # Models generated from the editor need to extract from the content field.
    json_data = json_data[:content]
  end
  load_solver(json_data)
end

"Parse and load a specific solver given the data as a dictionary."
function load_solver(json_data::Dict)::ReconciliationSolver
  load_solver(JSON3.read(JSON3.write(json_data)))
end

"Parse and load a specific solver given the data in JSON format."
function load_solver(json_data::JSON3.Object)::ReconciliationSolver
  solver_option = get(json_data, "solver", "")

  if solver_option == "Default"
    # Defaults to presolve + starting values + JuMP.
    DefaultSolver()

  elseif solver_option == "JuMP"
    JuMPSolver()

  else
    throw(ArgumentError("Solver $(solver_option) unknown."))
  end
end

function to_dict(obj::Any)
  if obj isa Dictionary
    Dict(k => to_dict(v) for (k, v) in pairs(obj))
  elseif obj isa Dict
    Dict(k => to_dict(v) for (k, v) in obj)
  elseif obj isa Vector
    [to_dict(entry) for entry in obj]
  else
    obj
  end
end


"Keys for additional serialization infos"
const additional_info_keys = [
  INFO_HEURISTIC,
  INFO_PRESOLVE_TIME,
  INFO_SV_TIME,
  INFO_SV_SOLVER,
  INFO_SV_SUCCEEDED,
  INFO_SV_PREP_TIME,
  INFO_SV_SOLVE_TIME,
  INFO_SV_NLS_TIME,
  INFO_SOLVE_TIME,
  INFO_JUMP_SOLVE_TIME,
  INFO_SOLVE2_TIME,
  INFO_JUMP_SOLVE2_TIME,
  INFO_JUMP_OPTIMIZER,
  INFO_SOLVER_TOLERANCE,
  INFO_MAX_ITERATIONS,
  INFO_JUMP_STATUS,
]

"Keys for serialization infos if solution is INFEASIBLE"
const relevant_error_keys = [
  INFO_ERROR_MESSAGE,
  INFO_REASON_INFEASIBILITY,
  INFO_PROBLEMATIC_REPLACEMENTS,
  INFO_CONFLICTING_REPLACEMENTS,
  INFO_CONFLICTING_EQ_DEPS,
  INFO_FAILURE_POINT,
  INFO_PREVIOUS_REPLACEMENTS,
]

"""
Extract solution data and convert to a dictionary.
Add more information from the solution to be serialized as needed.
"""
function to_dict(sol::ReconciliationSolution)::Dict
  # Make sure that the uncertainties field is present. If it was not calculated,
  # then it is stored as an empty field.
  uncertainties = if haskey(sol.infos, INFO_UNCERTAINTIES)
    uncertainties_dict = sol.infos[INFO_UNCERTAINTIES]
    to_dict(uncertainties_dict)
  else
    ""
  end

  # Extract information on unobservable variables if available
  removed_unobservable = if get(sol.infos, INFO_UNOBSERVABLE_REMOVED, false)
    "(removed)"
  else
    "(not removed)"
  end
  unobservable_vars_list = collect(keys(sol.unobservable_variables))
  unobservable_vars = if haskey(sol.infos, INFO_N_UNOBSERVABLE)
    string(sol.infos[INFO_N_UNOBSERVABLE], " unobservable unmeasured variables identified ", removed_unobservable, ": ")
  else
    "No unobservable variables identified or removed"
  end

  variables_value_dict = to_dict(sol.variables_value)
  equations_value_dict = to_dict(sol.score.equations_value)
  residual_value = if isnothing(sol.score.residual)
    nothing
  else
    round(sol.score.residual, sigdigits = 4)
  end
  ver = string(Pkg.project().version)

  returndict = Dict(
    "objective_value" => sol.objective_value,
    "variables_value" => variables_value_dict,
    "equations_value" => equations_value_dict,
    "residual_value" => isnothing(residual_value) ? "N/A" : residual_value,
    "uncertainties" => uncertainties,
    "status" => sol.status,
    "globality_proven" => sol.is_global,
    "version" => ver,
    "unobservable unmeasured variables" => [unobservable_vars, unobservable_vars_list],
  )

  # Serialize additional information if available
  merge!(returndict, Dict([k => sol.infos[k] for k in additional_info_keys if haskey(sol.infos, k)]))

  # Error information if INFEASIBLE
  if sol.status == INFEASIBLE || sol.status == CONTRADICTORY
    returndict["error information"] =
      Dict([k => to_dict(sol.infos[k]) for k in relevant_error_keys if haskey(sol.infos, k)])
  end

  return returndict
end

"Store the solution in JSON format."
function serialize(sol::ReconciliationSolution, file::String = tempname() * "sol.json")::String
  result = to_dict(sol)
  open(file, "w") do io
    JSON3.pretty(io, result)
  end
  file
end

# Helper to serialize bounds: converts missing to null and Inf to string representation
function serialize_bound(bound::Union{Float64, Missing})
  if ismissing(bound)
    nothing
  elseif bound == -Inf
    "-Inf"
  elseif bound == Inf
    "+Inf"
  else
    bound
  end
end

# Needs special casing because Inf and missing cannot be handled in JSON.
function to_dict(var::MeasuredVariable)::Dict
  Dict(
    "type" => "MeasuredVariable",
    "name" => var.name,
    "value" => var.value,
    "uncertainty" => var.uncertainty,
    "info" => var.info,
    "units" => ismissing(var.units) ? "" : var.units,
    "is_transfer_coefficient" => var.tc,
    "lb" => serialize_bound(var.lb),
    "ub" => serialize_bound(var.ub),
  )
end

function to_dict(var::UnmeasuredVariable)::Dict
  Dict(
    "type" => "UnmeasuredVariable",
    "name" => var.name,
    "info" => var.info,
    "units" => ismissing(var.units) ? "" : var.units,
    "is_transfer_coefficient" => var.tc,
    "lb" => serialize_bound(var.lb),
    "ub" => serialize_bound(var.ub),
  )
end
function to_dict(var::FixedVariable)::Dict
  Dict(
    "type" => "FixedVariable",
    "name" => var.name,
    "value" => var.value,
    "info" => var.info,
    "units" => ismissing(var.units) ? "" : var.units,
    "is_transfer_coefficient" => var.tc,
  )
end

function to_dict(rec::ReconciliationProblem)::Dict
  all_equations = [values(rec.equations)...]
  all_variables = [to_dict(get_variable(rec, name)) for name in keys(rec.variables)]
  Dict(
    "name" => rec.name,
    "equations" => all_equations,
    "variables" => all_variables,
  )
end

"Store the problem in JSON format."
function serialize(rec::ReconciliationProblem, file::String = tempname() * "problem.json")::String
  result = to_dict(rec)
  open(file, "w") do io
    JSON3.pretty(io, result)
  end
  file
end

end
