using Test
using OpenSEFA
using OpenSEFA.STANModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.JuMPInterfaceModule

# Filepath for each STAN model available under the `models` folder.
# Returns empty when MODELS_PATH doesn't exist; model-dependent tests skip gracefully.
models = all_models()

@testset "Solve models using Ipopt NLP via JuMP." begin
  for m in models
    @info "Solve using Ipopt NLP via JuMP - Running model : $(basename(m))"
    # Convert data to reconciliation problem.
    trace_filename = joinpath(m, "data", "trace.txt")
    if !isfile(trace_filename)
      @warn "Skipping model $(basename(m)). Expected trace file under `/data/trace.txt`."
      continue
    end
    # Convert data to reconciliation problem.
    strace = parse(STANTrace, m)
    rec = convert(ReconciliationProblem, strace)

    # Convert STAN results to reconciliation solution.
    sol = solve(rec, JuMPSolver())
    @test sol isa ReconciliationSolution
  end
end

@testset "Check presolve followed by JuMP solver." begin
  for m in models
    @info "Presolve followed by JuMP solver - Running model : $(basename(m))"
    # Convert data to reconciliation problem.
    strace = parse(STANTrace, m)
    rec = convert(ReconciliationProblem, strace)
    solver = DefaultSolver(; optimization_solver = JuMPSolver(), starting_values_solver = nothing)

    sol = solve(rec, solver)
    @test sol isa ReconciliationSolution
  end
end

# The solution can be computed analytically, see documentation.
# The system has two local minima which are also global at F1 = 1/5 and F1 = 4/5.
function uniqueness_1()
  # Define the equations based on the flowsheet example using string parsing.
  f1 = parse_equation("t13 + t14 = 1")
  f2 = parse_equation("t23 + t24 = 1")
  f3 = parse_equation("F3 = t13 * F1 + t23 * F2")
  f4 = parse_equation("F4 = t14 * F1 + t24 * F2")

  # Define the measured variables with their uncertainties.
  F1 = MeasuredVariable(; name = "F1", value = 1.0, uncertainty = 5 / 3)
  F2 = MeasuredVariable(; name = "F2", value = 4 / 25, uncertainty = 5 / 4)
  # This is a transfer coefficient.
  t13 = MeasuredVariable(; name = "t13", value = 1.0, uncertainty = 1.0, tc = true)

  # Define the unknown variables (not measured).
  t14 = UnmeasuredVariable(; name = "t14")
  t23 = UnmeasuredVariable(; name = "t23")

  # Define the variables which do not have uncertainty (fixed).
  F3 = FixedVariable(; name = "F3", value = 4 / 25)
  F4 = FixedVariable(; name = "F4", value = 1.0)
  t24 = FixedVariable(; name = "t24", value = 1.0)

  # Instantiate the reconciliation problem.
  rec = ReconciliationProblem(; name = "Non-uniqueness example")
  add_variable!(rec, [F1, F2, F3, F4, t13, t14, t23, t24])
  add_equation!(rec, [f1, f2, f3, f4])

  rec
end

# The solution can be computed analytically, see several_local_minima.pdf
# The system has two local minima at F1 = 14/9 and F1 = 22 but only the minimum at F1 = 22 is global.
function uniqueness_2()
  # Define the equations based on the flowsheet example using string parsing.
  f1 = parse_equation("t13 + t14 = 1")
  f2 = parse_equation("t23 + t24 = 1")
  f3 = parse_equation("F3 = t13 * F1 + t23 * F2")
  f4 = parse_equation("F4 = t14 * F1 + t24 * F2")

  # Define the measured variables with their uncertainties.
  F1 = MeasuredVariable(; name = "F1", value = 97 / 4, uncertainty = 14 * sqrt(5))
  F2 = MeasuredVariable(; name = "F2", value = 158 / 99, uncertainty = 7 * sqrt(15))
  # This is a transfer coefficient.
  t13 = MeasuredVariable(; name = "t13", value = 31777 / 41580, uncertainty = 1.0, tc = true)

  # Define the unknown variables (not measured).
  t14 = UnmeasuredVariable(; name = "t14", tc = true)
  t23 = UnmeasuredVariable(; name = "t23", tc = true)

  # Define the variables which do not have uncertainty (fixed).
  F3 = FixedVariable(; name = "F3", value = 1.0)
  F4 = FixedVariable(; name = "F4", value = 22.0)
  t24 = FixedVariable(; name = "t24", value = 1.0)

  # Instantiate the reconciliation problem.
  rec = ReconciliationProblem(; name = "Non-uniqueness example")
  add_variable!(rec, [F1, F2, F3, F4, t13, t14, t23, t24])
  add_equation!(rec, [f1, f2, f3, f4])

  rec
end

# The solution can be computed analytically, see several_local_minima.pdf
# The system has two local minima at F1 = 14/9 and F1 = 22 but only the minimum at F1 = 22 is global.
function uniqueness_3()
  # Parse the equations constructed by the SEFA editor using string parsing.
  f1 = parse_equation("F1 + F2 - F3 - F4 = 0")
  f2 = parse_equation("F3 - F1 * t13 - F2 * t23 = 0")
  f3 = parse_equation("F4 - F1 * t14 - F2 * t24 = 0")
  f4 = parse_equation("- t13 - t14 = -1.0")
  f5 = parse_equation("- t23 - t24 = -1.0")

  # Define the measured variables with their uncertainties.
  F1 = MeasuredVariable(; name = "F1", value = 97 / 4, uncertainty = 14 * sqrt(5))
  F2 = MeasuredVariable(; name = "F2", value = 158 / 99, uncertainty = 7 * sqrt(15))
  # This is a transfer coefficient.
  t13 = MeasuredVariable(; name = "t13", value = 31777 / 41580, uncertainty = 1.0, tc = true)

  # Define the unknown variables (not measured).
  t14 = UnmeasuredVariable(; name = "t14")
  t23 = UnmeasuredVariable(; name = "t23")

  # Define the variables which do not have uncertainty (fixed).
  F3 = FixedVariable(; name = "F3", value = 1.0)
  F4 = FixedVariable(; name = "F4", value = 22.0)
  t24 = FixedVariable(; name = "t24", value = 1.0)

  # Instantiate the reconciliation problem.
  rec = ReconciliationProblem(; name = "Non-uniqueness example")
  add_variable!(rec, [F1, F2, F3, F4, t13, t14, t23, t24])
  add_equation!(rec, [f1, f2, f3, f4, f5])

  rec
end

# The solution can be computed analytically, see several_local_minima.pdf
# The system has two local minima at F1 = 14/9 and F1 = 22 but only the minimum at F1 = 22 is global.
function uniqueness_4()
  # Load the model from a file.
  # The source file lives in MaterialFlowAnalysis/workbench; reference it via relative path.
  examples_dir = joinpath(pkgdir(OpenSEFA), "test", "examples")
  filepath = joinpath(examples_dir, "several_local_minima.json")
  @assert isfile(filepath)
  load_model(filepath)
end

# Check if the solution matches one of the expected analytic solutions.
function matches_solution(sol::ReconciliationSolution, expected::Dict{Label, Float64}; atol = 1e-5)
  vals = sol.variables_value
  all(var -> relaxed_isapprox(vals[var], expected[var], atol = atol), keys(expected))
end

@testset "Non-uniqueness example" begin
  # Create the reconciliation problem
  model = uniqueness_1()

  # Expected analytic solutions from non-uniqueness-example.pdf
  expected_sol1 = Dict("F1" => 1 / 5, "F2" => 24 / 25, "t13" => 4 / 5, "t14" => 1 / 5, "t23" => 0.0)
  expected_sol2 = Dict("F1" => 4 / 5, "F2" => 9 / 25, "t13" => 1 / 5, "t14" => 4 / 5, "t23" => 0.0)

  # The solvers find at least one of the optimal solutions.
  sol_default = solve(model)
  @test sol_default.status == FEASIBLE
  @test matches_solution(sol_default, expected_sol1) || matches_solution(sol_default, expected_sol2)

  sol_jump = solve(model, JuMPSolver())
  @test sol_jump.status == FEASIBLE
  @test matches_solution(sol_jump, expected_sol1) || matches_solution(sol_jump, expected_sol2)
end

@testset "Example with several local minima" begin
  @testset "using equations from mathematical model" begin
    # Create the reconciliation problem
    model = uniqueness_2()

    # Expected analytic solutions from several_local_minima.pdf
    expected_sol = Dict("F1" => 22.0, "F2" => 1.0, "t13" => 1 / 22, "t14" => 21 / 22, "t23" => 0.0)
    local_min = Dict("F1" => 14 / 9, "F2" => 193 / 9, "t13" => 9 / 14, "t14" => 5 / 14, "t23" => 0.0)

    # The solvers find at least one of the optimal solutions.
    sol_default = solve(model)
    @test sol_default.status == FEASIBLE
    @test matches_solution(sol_default, expected_sol)
    if matches_solution(sol_default, local_min)
      @warn "The DefaultSolver found the local (but not global) minimum as solution."
    end

    sol_jump = solve(model, JuMPSolver())
    @test sol_jump.status == FEASIBLE
    @test matches_solution(sol_jump, expected_sol)
    if matches_solution(sol_jump, local_min)
      @warn "The JuMPSolver found the local (but not global) minimum as solution."
    end
  end

  # Run again with the model as it is generated by the SEFA editor
  @testset "using equations as generated from SEFA Editor" begin
    # Create the reconciliation problem
    model = uniqueness_3()

    # Expected analytic solutions from several_local_minima.pdf
    expected_sol = Dict("F1" => 22.0, "F2" => 1.0, "t13" => 1 / 22, "t14" => 21 / 22, "t23" => 0.0)
    local_min = Dict("F1" => 14 / 9, "F2" => 193 / 9, "t13" => 9 / 14, "t14" => 5 / 14, "t23" => 0.0)

    # The solvers find at least one of the optimal solutions.
    sol_default = solve(model)
    @test sol_default.status == FEASIBLE broken = false
    @test matches_solution(sol_default, expected_sol) broken = true
    if matches_solution(sol_default, local_min)
      @warn "The DefaultSolver found the local (but not global) minimum as solution."
    end

    sol_jump = solve(model, JuMPSolver())
    @test sol_jump.status == FEASIBLE
    @test matches_solution(sol_jump, expected_sol)
    if matches_solution(sol_jump, local_min)
      @warn "The JuMPSolver found the local (but not global) minimum as solution."
    end
  end

  # Run again with the model generated by the SEFA editor
  @testset "using model from SEFA Editor" begin
    # Create the reconciliation problem
    model = uniqueness_4()

    # Expected analytic solutions from several_local_minima.pdf
    expected_sol = Dict(
      "F1_mass_N1" => 22.0,
      "F2_mass_N1" => 1.0,
      "P1_mass_F1F3_N1" => 1 / 22,
      "P1_mass_F1F4_N1" => 21 / 22,
      "P1_mass_F2F3_N1" => 0.0,
    )
    local_min = Dict(
      "F1_mass_N1" => 14 / 9,
      "F2_mass_N1" => 193 / 9,
      "P1_mass_F1F3_N1" => 9 / 14,
      "P1_mass_F1F4_N1" => 5 / 14,
      "P1_mass_F2F3_N1" => 0.0,
    )

    # The solvers find at least one of the optimal solutions.
    sol_default = solve(model)
    @test sol_default.status == FEASIBLE
    @test matches_solution(sol_default, expected_sol) broken = true
    if matches_solution(sol_default, local_min)
      @warn "The DefaultSolver found the local (but not global) minimum as solution."
    end

    sol_jump = solve(model, JuMPSolver())
    @test sol_jump.status == FEASIBLE
    @test matches_solution(sol_jump, expected_sol)
    if matches_solution(sol_jump, local_min)
      @warn "The JuMPSolver found the local (but not global) minimum as solution."
    end
  end
end
