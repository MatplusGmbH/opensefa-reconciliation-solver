using Test
import OpenSEFA
using OpenSEFA.NonlinearSolveInterfaceModule
using OpenSEFA.ConstraintsModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.DefaultSolverModule

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))

@testset "NonlinearSolver constructor" begin
  solver = NonlinearSolver()
  @test solver isa NonlinearSolver
end

@testset "NonlinearSolver solves seven flows" begin
  rec = seven_flows_cencic_2012()
  sol = solve(rec, NonlinearSolver())
  @test sol isa ReconciliationSolution
  @test sol.score.residual < 1e-4
end

@testset "NonlinearSolver pipeline init and solve!" begin
  rec = seven_flows_cencic_2012()
  pipe = init(rec, NonlinearSolver())
  sol = solve!(pipe)
  @test sol isa ReconciliationSolution
end
