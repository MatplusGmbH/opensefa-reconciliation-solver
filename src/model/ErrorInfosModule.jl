"""
Infos infrastructure for `ReconciliationSolution.infos`.

This module owns:
- All `INFO_*` string key constants — the single authoritative source for every key
  written into a `Dictionary{String, Any}` infos payload.
- Helper functions for reading and writing infos dictionaries (`add_info!`,
  `merge_info!`, `overwrite_info!`, `error_infos`).

Placing these here (rather than in `ReconciliationModule`) keeps the heavyweight
algorithms module free of boilerplate and makes the constants available to the IO
layer (`STANModule`) which is loaded before `ReconciliationModule`.
"""
module ErrorInfosModule

using Dictionaries
using ..ErrorModule

export
  # Helpers
  add_info!,
  merge_info!,
  overwrite_info!,
  error_infos,
  # Info key constants
  INFO_PRESOLVE_TIME,
  INFO_SV_TIME,
  INFO_SV_PREP_TIME,
  INFO_SV_SOLVE_TIME,
  INFO_SV_NLS_TIME,
  INFO_SOLVE_TIME,
  INFO_SOLVE2_TIME,
  INFO_JUMP_SOLVE_TIME,
  INFO_JUMP_SOLVE2_TIME,
  INFO_STAN_MODELLING_TIME,
  INFO_STAN_COMPUTATION_TIME,
  INFO_HEURISTIC,
  INFO_SV_SOLVER,
  INFO_SV_SUCCEEDED,
  INFO_JUMP_OPTIMIZER,
  INFO_JUMP_STATUS,
  INFO_SOLVER_TOLERANCE,
  INFO_MAX_ITERATIONS,
  INFO_NLS_RETCODE,
  INFO_NLS_RESID,
  INFO_REDUCED_PROBLEM,
  INFO_JUMP_MODEL,
  INFO_COMPUTED_VARIABLES,
  INFO_MISSING_VARIABLES,
  INFO_COMPUTED_EQUATIONS,
  INFO_MISSING_EQUATIONS,
  INFO_COMPUTED_VARIABLES_EXTREMA,
  INFO_TC_EXTREMA,
  INFO_N_UNOBSERVABLE,
  INFO_UNOBSERVABLE_REMOVED,
  INFO_UNCERTAINTIES,
  INFO_DUAL,
  INFO_REASON_INFEASIBILITY,
  INFO_FAILURE_POINT,
  INFO_ERROR_MESSAGE,
  INFO_PROBLEMATIC_REPLACEMENTS,
  INFO_CONFLICTING_REPLACEMENTS,
  INFO_CONFLICTING_EQ_DEPS,
  INFO_PREVIOUS_REPLACEMENTS,
  INFO_STAN_STATUS,
  INFO_STAN_VERSION,
  INFO_STAN_MACHINE

# ====================
# Info key constants
# ====================
# All keys written into ReconciliationSolution.infos must be declared here.
# Using constants instead of raw string literals makes typos a grep-time error
# and enables safe cross-file renames.

# Timing
const INFO_PRESOLVE_TIME              = "Presolve time (s)"
const INFO_SV_TIME                    = "Starting values time (s)"
const INFO_SV_PREP_TIME               = "Starting values preparation time (s)"
const INFO_SV_SOLVE_TIME              = "Starting values solving time (s)"
const INFO_SV_NLS_TIME                = "Starting values NonlinearSolve time (s)"
const INFO_SOLVE_TIME                 = "Solve time (s)"
const INFO_SOLVE2_TIME                = "Solve second time (s)"
const INFO_JUMP_SOLVE_TIME            = "JuMP solve time (s)"
const INFO_JUMP_SOLVE2_TIME           = "JuMP second solve time (s)"
const INFO_STAN_MODELLING_TIME        = "Modelling time (s)"
const INFO_STAN_COMPUTATION_TIME      = "Computation time (s)"

# Solver identity / control
const INFO_HEURISTIC                  = "Heuristic"
const INFO_SV_SOLVER                  = "Starting values solver"
const INFO_SV_SUCCEEDED               = "Starting values succeeded"
const INFO_JUMP_OPTIMIZER             = "JuMP optimizer"
const INFO_JUMP_STATUS                = "JuMP status"
const INFO_SOLVER_TOLERANCE           = "Solver tolerance"
const INFO_MAX_ITERATIONS             = "Max. Iterations"

# NLS retcode / residual (NonlinearSolveInterfaceModule)
const INFO_NLS_RETCODE                = "retcode"
const INFO_NLS_RESID                  = "resid"

# Problem structure
const INFO_REDUCED_PROBLEM            = "Reduced problem"
const INFO_JUMP_MODEL                 = "Model"
const INFO_COMPUTED_VARIABLES         = "Computed variables"
const INFO_MISSING_VARIABLES          = "Missing variables"
const INFO_COMPUTED_EQUATIONS         = "Computed equations"
const INFO_MISSING_EQUATIONS          = "Missing equations"
const INFO_COMPUTED_VARIABLES_EXTREMA = "Computed variables extrema"
const INFO_TC_EXTREMA                 = "Transfer coefficients extrema"

# Unobservability
const INFO_N_UNOBSERVABLE             = "Number of unobservable variables"
const INFO_UNOBSERVABLE_REMOVED       = "Unobservable variables removed"

# Error propagation
const INFO_UNCERTAINTIES              = "Uncertainties"
const INFO_DUAL                       = "Dual"

# Failure diagnostics
const INFO_REASON_INFEASIBILITY       = "Reason for infeasibility"
const INFO_FAILURE_POINT              = "Failure point"
const INFO_ERROR_MESSAGE              = "Error message"
const INFO_PROBLEMATIC_REPLACEMENTS   = "Problematic replacements"
const INFO_CONFLICTING_REPLACEMENTS   = "Conflicting replacements"
const INFO_CONFLICTING_EQ_DEPS        = "Conflicting equation dependencies"
const INFO_PREVIOUS_REPLACEMENTS      = "Previous replacements"

# STAN integration
const INFO_STAN_STATUS                = "STAN status"
const INFO_STAN_VERSION               = "STAN version"
const INFO_STAN_MACHINE               = "Machine information"

# =================
# Helper functions
# =================

"Add an info to a dictionary. Check if the info key already exists."
function add_info!(dict::Dictionary{String, Any}, key::String, new_info::Any)
  value_old = get(dict, key, nothing)
  if isnothing(value_old)
    value_new = new_info
  elseif value_old isa Vector{String}
    value_new = push!(value_old, new_info)
  else
    value_new = [value_old, new_info]
  end
  set!(dict, key, value_new)
end

"Merge infos from two dictionaries. Optionally provide a list of keys that shall be merged from dict2."
function merge_info!(
  dict1::Dictionary{String, Any},
  dict2::Dictionary{String, Any},
  keylist::Union{Vector{String}, Nothing} = nothing,
)
  if isnothing(keylist)
    keylist = keys(dict2)
  else
    filter!(el -> el in keys(dict2), keylist)
  end
  for key in keylist
    add_info!(dict1, key, dict2[key])
  end
end

"""
Overwrite entries in `dict1` with values from `dict2` for the given keys.

Unlike `merge_info!`, which accumulates multiple values into a vector when a key is
already present, this function always replaces the existing value with `set!`. Use it
for structural results (e.g. `INFO_REDUCED_PROBLEM`, `INFO_UNCERTAINTIES`, `INFO_DUAL`)
that are single-valued by design and must not accumulate across solver stages or retry
attempts.
"""
function overwrite_info!(
  dict1::Dictionary{String, Any},
  dict2::Dictionary{String, Any},
  keylist::Vector{String},
)
  for key in keylist
    if haskey(dict2, key)
      set!(dict1, key, dict2[key])
    end
  end
end

"""
Return the infos dictionary attached to a `SolverError`, if any.

The default implementation returns an empty dictionary, so callers never need to
check for the field's existence with `hasfield`. Override for error types that
carry an `infos` field.
"""
error_infos(::SolverError) = Dictionary{String, Any}()
error_infos(e::InfeasibleEquationError) = e.infos
error_infos(e::IllposedProblemError) = e.infos

end
