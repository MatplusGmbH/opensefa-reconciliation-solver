using Test
using ResultTypes
import OpenSEFA
using OpenSEFA.DefaultSolverModule
using OpenSEFA.JuMPInterfaceModule
using OpenSEFA.ConstraintsModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.PresolveModule
using OpenSEFA.ErrorModule
using OpenSEFA.SerializationModule

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))

@testset "DefaultSolver constructor" begin
  solver = DefaultSolver()
  @test solver isa DefaultSolver
end

@testset "DefaultSolver solves seven flows" begin
  rec = seven_flows_cencic_2012()
  sol = solve(rec, DefaultSolver())
  @test sol isa ReconciliationSolution
  @test sol.status in SUCCESS_CODES
  @test sol.score.residual < 1e-6
end

@testset "DefaultSolver with custom optimization_solver" begin
  rec = seven_flows_cencic_2012()
  sol = solve(rec, DefaultSolver(; optimization_solver = JuMPSolver()))
  @test sol isa ReconciliationSolution
  @test sol.status in SUCCESS_CODES
end

@testset "DefaultPipeline init and solve!" begin
  rec = seven_flows_cencic_2012()
  pipe = init(rec, DefaultSolver())
  @test pipe isa DefaultPipeline
  sol = solve!(pipe)
  @test sol isa ReconciliationSolution
  @test sol.status in SUCCESS_CODES
end

@testset "find_starting_values" begin
  rec = seven_flows_cencic_2012()
  remove_fixed_variables!(rec)
  result = find_starting_values(rec)
  # Problem is feasible with unobservable variables.
  @test result isa Result{StartingReconciliationProblem, SolverError}
end

@testset "Starting value solution only used when basic constraints are satisfied" begin
  # Regression test: when the starting-value solution violates basic constraints (e.g.
  # lb/ub bounds), DefaultSolver must discard it and proceed to the optimization stage
  # rather than returning a solution that breaks the bounds.
  filepath = joinpath(pkgdir(OpenSEFA), "test", "examples", "feasability_issue.json")
  @assert isfile(filepath)
  model = load_model(filepath)

  sol1 = solve(model; force_basic_constraints = false) # without enforcement
  sol2 = solve(model)                                  # with enforcement (default)
  if has_basic_constraints(sol1)
    @warn "Meaningless test: basic constraints satisfied even without enforcement."
  end
  @test has_basic_constraints(sol2)
  @test sol2.status == FEASIBLE_WITH_UNOBSERVABLE_VARIABLES
end
