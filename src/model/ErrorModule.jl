"""
Custom exception types for the solver.

The error hierarchy is:
- `SolverError` — abstract base for all solver errors.
  - `ContradictionError` — the problem is contradictory (no feasible solution exists).
    - `InfeasibleEquationError` — a specific equation cannot be satisfied.
    - `IllposedProblemError` — the problem is ill-posed (e.g. duplicate variable names).
  - `InternalError` — unexpected internal state; indicates a bug.
- `FileIOError` — file input/output failure (not a `SolverError`).
"""
module ErrorModule

using Dictionaries

export
  SolverError,
  ContradictionError,
  InfeasibleEquationError,
  InternalError,
  IllposedProblemError,
  FileIOError

"Generic solver error."
abstract type SolverError <: Exception end

"Solver error due to a recognized contradiction."
abstract type ContradictionError <: SolverError end

"Error representing an infeasible equation in a set of equations."
struct InfeasibleEquationError <: ContradictionError
  "Error description."
  msg::String
  "Index of the infeasible equation."
  index::Vector{Int64}
  "Residual value."
  residual::Vector{Float64}
  "Additional information for debugging."
  infos::Dictionary{String, Any}
end

function InfeasibleEquationError(msg::String;
    index::Vector{Int64}=Int64[],
    residual::Vector{Float64}=Float64[],
    infos::Dictionary{String, Any}=Dictionary{String, Any}())
  InfeasibleEquationError(msg, index, residual, infos)
end

"Generic solver error."
struct InternalError <: SolverError
  "Error description."
  msg::String
end

"Error representing an illposed problem."
struct IllposedProblemError <: ContradictionError
  "Error description."
  msg::String
  "Name of the illposed problem."
  name::String
  "Additional information for debugging."
  infos::Dictionary{String, Any}
end

function IllposedProblemError(msg::String;
    name::String="",
    infos::Dictionary{String, Any}=Dictionary{String, Any}())
  IllposedProblemError(msg, name, infos)
end

"File input/output error."
struct FileIOError <: Exception
  "Error description."
  msg::String
end

end
