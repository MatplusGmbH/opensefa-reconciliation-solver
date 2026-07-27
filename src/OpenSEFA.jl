# Copyright (c) Matplus GmbH
# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0

"""
OpenSEFA is a Julia package for Material Flow Analysis (MFA) and
Substance & Energy Flow Analysis (SEFA) via data reconciliation,
distributed under the PolyForm Noncommercial License 1.0.0.

It provides:
- A constraint-based model for material/substance flow balances.
- Data reconciliation algorithms to find statistically optimal values consistent with model constraints.
- A presolver that simplifies problems through structural reductions before optimization.
- Multiple solver backends: Ipopt, Alpine, and NonlinearSolve.
- I/O support for STAN, JSON, and Excel formats.
- An HTTP server exposing a REST API for remote solving.
"""
module OpenSEFA

using DotEnv
using PrecompileTools
using Reexport

function __init__()
  path = joinpath(dirname(@__DIR__), ".env")
  if isfile(path)
    DotEnv.load!(path)
  end
end

include("sources.jl")

# External re-exports
@reexport using Dictionaries
@reexport using CommonSolve: solve
@reexport using LinearAlgebra
@reexport using ResultTypes

# Internal modules
@reexport using .ErrorModule
@reexport using .ErrorInfosModule
@reexport using .ComparisonModule
@reexport using .ConstraintsModule
@reexport using .ConstraintParserModule
@reexport using .RowEchelonModule
@reexport using .SparseRowEchelonModule
@reexport using .ReconciliationModule
@reexport using .PresolveModule
@reexport using .QCQPModule
@reexport using .JuMPInterfaceModule
@reexport using .IpoptSolverModule
@reexport using .NonlinearSolveInterfaceModule
@reexport using .DefaultSolverModule
@reexport using .STANModule
@reexport using .SerializationModule
@reexport using .ExcelModule
@reexport using .HTTPServerModule

include("precompilation.jl")

end # module OpenSEFA
