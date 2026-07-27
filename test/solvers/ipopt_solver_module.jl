using Test
import OpenSEFA 
using OpenSEFA.IpoptSolverModule
using OpenSEFA.JuMPInterfaceModule
using OpenSEFA.ConstraintsModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.DefaultSolverModule

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))

@testset "IpoptSolver constructor" begin
  solver = IpoptSolver()
  @test solver isa JuMPSolver
end

@testset "IpoptSolver accepts keyword arguments" begin
  solver = IpoptSolver(; max_iter = 100)
  @test solver isa JuMPSolver
  @test solver.max_iter == 100
end

@testset "IpoptSolver solves seven flows" begin
  rec = seven_flows_cencic_2012()
  sol = solve(rec, IpoptSolver())
  @test sol isa ReconciliationSolution
  @test sol.status in SUCCESS_CODES
  @test sol.score.residual < 1e-6
end
