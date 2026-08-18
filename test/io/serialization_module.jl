using Test
using JSON3
using Dictionaries
import OpenSEFA
using OpenSEFA.SerializationModule
using OpenSEFA.DefaultSolverModule
using OpenSEFA.JuMPInterfaceModule
using OpenSEFA.NonlinearSolveInterfaceModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.STANModule
using OpenSEFA.ComparisonModule
using OpenSEFA.ConstraintsModule

EXAMPLES_DIR = joinpath(pkgdir(OpenSEFA), "test", "examples")

@testset "Load model defined as JSON." begin
  example_json = joinpath(EXAMPLES_DIR, "seven_flows.json")
  rec = load_model(example_json)
  @test length(rec.equations) == 4
  @test length(rec.variables) == 8
end

@testset "Load variables and equations from a JSON model." begin
  example_json = joinpath(EXAMPLES_DIR, "serialization_test_model.json")
  rec = load_model(example_json)

  vars_expected =
    dictionary(["v" => MeasuredVariable, "w" => FixedVariable, "x" => MeasuredVariable, "y" => UnmeasuredVariable])
  @test rec.variables == vars_expected
  @test rec.equations[1] ==
        Equation(ConstantTerm(0.0), LinearTerm[LinearTerm("w", 1.0)], BilinearTerm[BilinearTerm("v", "x", -1.0)])
  @test rec.equations[2] == Equation(ConstantTerm(-9.0), LinearTerm[], BilinearTerm[BilinearTerm("v", "x", 3.0)])
  @test rec.equations[3] == Equation(ConstantTerm(6.0), LinearTerm[LinearTerm("w", -2.0)], BilinearTerm[])
end

@testset "Store solution." begin
  example_json = joinpath(EXAMPLES_DIR, "seven_flows.json")
  rec = load_model(example_json)
  sol = solve(rec)
  file = serialize(sol)

  # Check conversion.
  json_sol = JSON3.read(file)
  expected = JSON3.read(
    JSON3.write(
      Dict(
        "variables_value" => Dict(
          "tc34" => 0.49597295051198165,
          "m6" => 76.21297998504508,
          "m3" => 302.416229692811,
          "m2" => 50,
          "m4" => 149.99026972340005,
          "m5" => 152.4259599700902,
          "m1" => 102.42595947011199,
          "m7" => 76.21297998504511,
        ),
        "objective_value" => 0.29591275484549834,
        "residual_value" => 4.993e-7,
      ),
    ),
  )

  # Test the objective value
  @test json_sol.objective_value ≈ expected.objective_value atol = 5e-7

  # Test the residual value
  @test json_sol.residual_value ≈ expected.residual_value atol = 5e-7

  # Test the variables' values
  for (var, expected_val) in expected.variables_value
    @test relaxed_isapprox(json_sol.variables_value[var], expected_val, atol = 5e-7)
  end
end

@testset "Store a model." begin
  # Check the roundtrip: loading a model from JSON, serializing it and then loading again
  # should match the original model.
  example_json = joinpath(EXAMPLES_DIR, "seven_flows.json")
  rec = load_model(example_json)
  aux_json = serialize(rec)
  rec_loaded = load_model(aux_json)
  @test rec_loaded.variables == rec.variables
  @test rec_loaded.equations == rec.equations
end

@testset "Different solver options." begin
  json_path = joinpath(EXAMPLES_DIR, "simple_two_flows_auto.json")
  model = load_model(json_path)
  solver = load_solver(json_path)
  @test typeof(solver) == DefaultSolver{NonlinearSolver, JuMPSolver, OpenSEFA.PresolveModule.MixConfig}
  sol = solve(model, solver)
  @test sol.status == FEASIBLE

  json_path = joinpath(EXAMPLES_DIR, "simple_two_flows_JuMP.json")
  model = load_model(json_path)
  solver = load_solver(json_path)
  @test typeof(solver) == JuMPSolver
  sol = solve(model, solver)
  @test sol.status == FEASIBLE

  json_path = joinpath(EXAMPLES_DIR, "simple_two_flows_JuMPpre.json")
  model = load_model(json_path)
  solver = load_solver(json_path)
  @test typeof(solver) == JuMPSolver
  sol = solve(model, solver)
  @test sol.status == FEASIBLE
end

@testset "Logging error information for infeasible problems" begin
  if isdir(MODELS_PATH)
    m = only(filter(x -> basename(x) == "Contradiction1_STAN", all_models()))
    strace = parse(STANTrace, m)
    rec = convert(ReconciliationProblem, strace)
    sol = solve(rec)
    # Check that the returned dict contains error information.
    returndict = to_dict(sol)
    @test haskey(returndict, "error information")
  else
    @warn "Skipping infeasible model test: MODELS_PATH not found."
  end
end

@testset "Converting nested Dictionary objects to Dict for serialization" begin
  # A Dictionary{String,Float64} nested as a value inside another Dictionary —
  # this mirrors how "Previous replacements" (Dictionary{Label,Float64}) is
  # stored inside sol.infos (Dictionary{String,Any}).
  inner = Dictionary(["x", "y"], [1.0, 2.0])
  outer = Dictionary(["nested", "scalar"], Any[inner, 42])

  result = to_dict(outer)

  # Top level must be a plain Dict, not a Dictionary.
  @test result isa Dict

  # The nested value must also have been converted.
  @test result["nested"] isa Dict
  @test result["nested"]["x"] == 1.0
  @test result["nested"]["y"] == 2.0

  # Non-Dictionary values must be left unchanged.
  @test result["scalar"] == 42

  # Same check for a Vector wrapping a Dictionary (mirrors "Failure point"
  # which can be a Vector{String}, but also covers Vector{Dictionary} paths).
  vec_outer = Dictionary(["list"], Any[[inner, inner]])
  vec_result = to_dict(vec_outer)
  @test vec_result isa Dict
  @test vec_result["list"] isa Vector
  @test all(d -> d isa Dict, vec_result["list"])
end
