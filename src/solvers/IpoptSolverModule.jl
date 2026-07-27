"""
Convenience module exposing `IpoptSolver()` for data reconciliation problems.

`IpoptSolver` is a named constructor for `JuMPSolver` that makes the optimizer choice
explicit. It is suitable for smooth (possibly nonlinear) reconciliation problems. For
nonconvex problems consider `AlpineSolver()` from `AlpineSolverModule`.
"""
module IpoptSolverModule

import Ipopt

using ..JuMPInterfaceModule

export IpoptSolver

"""
Convenience constructor for a `JuMPSolver` using `Ipopt.Optimizer`.
Accepts the same keyword arguments as `JuMPSolver`.
"""
IpoptSolver(; kwargs...) = JuMPSolver(; optimizer = Ipopt.Optimizer, kwargs...)

end
