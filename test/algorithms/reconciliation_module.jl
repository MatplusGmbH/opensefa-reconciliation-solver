using Test
using Suppressor
using Dictionaries
using LinearAlgebra
using SparseArrays
import OpenSEFA
using OpenSEFA.STANModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.ConstraintsModule
using OpenSEFA.PresolveModule
using OpenSEFA.NonlinearSolveInterfaceModule
using OpenSEFA.QCQPModule
using OpenSEFA.DefaultSolverModule
using OpenSEFA.RowEchelonModule

# Filepath for each STAN model available under the `models` folder.
models = all_models()

@testset "Adding and removing variables." begin
  rec = ReconciliationProblem()
  z1 = FixedVariable(; name = "z1", value = 1.0)

  add_variable!(rec, z1)
  @test isempty(rec.var_to_eq["z1"])
  @test rec.variables["z1"] == FixedVariable
  @test rec.fixed["z1"] == FixedVariable("z1", 1.0, "", "", false)

  remove_variable!(rec, "z1")
  @test isempty(rec.variables)
  @test isempty(rec.fixed)
  @test isempty(rec.var_to_eq)
end

@testset "Adding and removing equations." begin
  rec = ReconciliationProblem()
  z1 = FixedVariable(; name = "z1", value = 1.0)
  m1 = MeasuredVariable(; name = "m1", value = 10.0, uncertainty = 1.0)
  add_variable!(rec, z1)
  add_variable!(rec, m1)
  eq = Equation(LinearTerm[LinearTerm("z1", -1.0), LinearTerm("m1", 1.0)], BilinearTerm[])
  add_equation!(rec, eq)

  @test rec.equations[1] == eq
  @test rec.var_to_eq["z1"] == [1]
  @test rec.var_to_eq["m1"] == [1]
  @test rec.eq_to_var[1] == ["z1", "m1"]

  remove_equation!(rec, 1)
  @test isempty(rec.equations)
  # Removing the equation does not remove variables.
  @test !isempty(rec.variables)
  @test isempty(rec.eq_to_var)
  @test rec.var_to_eq["z1"] == EID[]
  @test rec.var_to_eq["m1"] == EID[]
end

@testset "Creating problems." begin
  @testset "Adding equation with all variables in the problem." begin
    z1 = FixedVariable(name = "z1", value = 2.0)
    z2 = FixedVariable(name = "z2", value = 3.0)
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 1.0)

    # Equation in the form of z1 - x2 = 0
    eq1 = Equation(LinearTerm[LinearTerm("z1", 1.0), LinearTerm("x2", -1.0)],
      BilinearTerm[],
    )

    # Equation in the form of z2 - x1 = 0
    eq2 = Equation(LinearTerm[LinearTerm("z2", 1.0), LinearTerm("x1", -1.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Testing add_equation")
    add_variable!(rec, [z1, z2, x1, x2])
    add_equation!(rec, [eq1, eq2])

    # Test if rec.equations contains eq1 and eq2
    @test rec.equations[1] == eq1
    @test rec.equations[2] == eq2
    @test rec.eq_to_var[1] == ["z1", "x2"]
    @test rec.eq_to_var[2] == ["z2", "x1"]
    @test rec.var_to_eq["z1"] == [1]
    @test rec.var_to_eq["z2"] == [2]
    @test rec.var_to_eq["x1"] == [2]
    @test rec.var_to_eq["x2"] == [1]
  end

  @testset "Adding equation with some variables not in the problem." begin
    z1 = FixedVariable(name = "z1", value = 2.0)
    z2 = FixedVariable(name = "z2", value = 3.0)
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 1.0)

    # Equation in the form of z1 - x2 = 0
    eq1 = Equation(LinearTerm[LinearTerm("z1", 1.0), LinearTerm("x2", -1.0)],
      BilinearTerm[],
    )

    # Equation in the form of z2 - x1 = 0
    eq2 = Equation(LinearTerm[LinearTerm("z2", 1.0), LinearTerm("x1", -1.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Testing add_equation")
    add_variable!(rec, z1)
    add_variable!(rec, z2)

    # ArgumentError thrown when adding an equation with variables not in rec.
    @test_throws ArgumentError add_equation!(rec, eq1)
    @test_throws ArgumentError add_equation!(rec, eq2)
  end

  @testset "Duplicate equations." begin
    z1 = FixedVariable(name = "z1", value = 2.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 1.0)

    # Equation in the form of z1 - x2 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("z1", 1.0), LinearTerm("x2", -1.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Testing duplicate equations")
    add_variable!(rec, z1)
    add_variable!(rec, x2)
    add_equation!(rec, eq1)
    # Add the same equation again.
    add_equation!(rec, eq1)
    @test length(rec.equations) == 2
    # Test that redundant equation was removed.
    remove_redundant_equations!(rec)
    @test length(rec.equations) == 1
    @test only(rec.equations) == eq1
  end

end

@testset "Unsafe add_equation (strict=false)." begin
  @testset "Adding equation with undefined variables." begin
    rec = ReconciliationProblem(name = "Testing undefined variables")
    eq1 = Equation(LinearTerm[LinearTerm("x1", 1.0)], BilinearTerm[])
    @test_throws ArgumentError add_equation!(rec, eq1)

    # Adding an equation with an undefined variable should work with strict=false
    add_equation!(rec, eq1, strict = false)
    @test length(rec.equations) == 1

    # Add the variables first, then the equation.
    rec = ReconciliationProblem(name = "Testing undefined variables")
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    add_variable!(rec, x1)
    add_equation!(rec, eq1)
    @test length(rec.equations) == 1
    @test length(rec.variables) == 1
  end

  @testset "Adding equation with repeated linear terms." begin
    rec = ReconciliationProblem(name = "Testing repeated linear terms")
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    add_variable!(rec, x1)

    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x1", 2.0)],
      BilinearTerm[],
    )

    @test_throws ArgumentError add_equation!(rec, eq1)

    # Should allow adding despite repeated linear terms.
    add_equation!(rec, eq1, strict = false)
    @test length(rec.equations) == 1
  end

  @testset "Adding equation with repeated bilinear terms" begin
    rec = ReconciliationProblem(name = "Testing repeated bilinear terms")
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 1.0)
    add_variable!(rec, x1)
    add_variable!(rec, x2)

    eq1 = Equation(
      LinearTerm[],
      BilinearTerm[BilinearTerm("x1", "x2", 1.0), BilinearTerm("x1", "x2", 2.0)],
    )

    @test_throws ArgumentError add_equation!(rec, eq1)

    # Should allow adding despite repeated bilinear terms
    add_equation!(rec, eq1, strict = false)
    @test length(rec.equations) == 1
  end

  @testset "Adding multiple invalid equations." begin
    rec = ReconciliationProblem(name = "Testing multiple invalid equations")
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    add_variable!(rec, x1)

    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x1", 2.0)],
      BilinearTerm[BilinearTerm("x1", "x1", 3.0), BilinearTerm("x1", "x1", 4.0)],
    )

    eq2 = Equation(
      LinearTerm[LinearTerm("y1", 1.0)],
      BilinearTerm[BilinearTerm("y1", "y2", 2.0)],
    )

    # Both should be allowed despite containing issues.
    add_equation!(rec, eq1, strict = false)
    add_equation!(rec, eq2, strict = false)
    @test length(rec.equations) == 2
  end
end

@testset "Converting the trace to reconciliation solution." begin
  for m in models
    # Convert data to reconciliation problem.
    trace_filename = joinpath(m, "data", "trace.txt")
    if !isfile(trace_filename)
      @warn "Skipping model $(basename(m)). Expected trace file under `/data/trace.txt`."
      continue
    end
    # Convert data to reconciliation problem.
    strace = parse(STANTrace, m)
    rec = convert(ReconciliationProblem, strace)
    @test rec isa ReconciliationProblem

    # Convert STAN results to reconciliation solution.
    @suppress_err begin
      sol_stan = convert(ReconciliationSolution, strace)
      @test sol_stan isa ReconciliationSolution
    end

    # Check results using seed array values.
    @suppress_err begin
      sol_stan_initial = convert(ReconciliationSolution, strace; use_initial_results = true)
      @test sol_stan_initial isa ReconciliationSolution
    end
  end
end

@testset "Solution score." begin
  @testset "Solution score for duplicate equations." begin
    problem = ReconciliationProblem(; name = "Problem with duplicate equation")
    eq = Equation(ConstantTerm(), LinearTerm[], BilinearTerm[])
    add_equation!(problem, eq)
    add_equation!(problem, eq)
    var = UnmeasuredVariable(; name = "var")
    add_variable!(problem, var)
    @test length(problem.equations) == 2
    sol = ReconciliationSolution(; problem, variables_value = Dictionary{Label, Float64}(["var"], [2.0]))
    @test sol.score.solved_variables_score == 1.0
    @test sol.score.solved_equations_score == 1.0
  end
  @testset "Solution score without equations." begin
    problem = ReconciliationProblem(; name = "Problem without equation")
    var = UnmeasuredVariable(; name = "var")
    add_variable!(problem, var)
    @test length(problem.equations) == 0
    sol = ReconciliationSolution(; problem, variables_value = Dictionary{Label, Float64}(["var"], [2.0]))
    @test sol.score.solved_variables_score == 1.0
    @test sol.score.solved_equations_score == 1.0
  end
end

@testset "Compare solutions" begin
  problem1 = ReconciliationProblem(; name = "Test problem 1")
  problem2 = ReconciliationProblem(; name = "Test problem 2")
  problem3 = ReconciliationProblem(; name = "Test problem 3")
  var1 = UnmeasuredVariable(; name = "var1")
  var2 = UnmeasuredVariable(; name = "var2")
  eq1 = Equation(ConstantTerm(-1.5), LinearTerm[LinearTerm("var1", 1.0)], BilinearTerm[])
  add_variable!(problem1, var1)
  add_variable!(problem2, var2)
  add_variable!(problem3, var1)
  add_equation!(problem3, eq1)
  sol1 = ReconciliationSolution(; problem = problem1, variables_value = Dictionary{Label, Float64}(["var1"], [1.0]))
  sol2 = ReconciliationSolution(; problem = problem2, variables_value = Dictionary{Label, Float64}(["var2"], [2.0]))
  sol3 = ReconciliationSolution(; problem = problem1, variables_value = Dictionary{Label, Float64}(["var1"], [3.0]))
  sol4 = ReconciliationSolution(; problem = problem3, variables_value = Dictionary{Label, Float64}(["var1"], [1.0]))
  sol5 = ReconciliationSolution(; problem = problem3, variables_value = Dictionary{Label, Float64}(["var1"], [3.0]))

  @testset "Compare solutions 1 and 2" begin
    # Since `sol1` and `sol2` hold different problems, it is expected to fail.
    @test_throws AssertionError compare(sol1, sol2)
  end

  @testset "Compare solutions 1 and 3" begin
    comp13 = compare(sol1, sol3)
    @test comp13.vars_absolute_deviation == Dictionary{Label, Float64}(["var1"], [2.0])
    @test comp13.vars_relative_deviation == Dictionary{Label, Float64}(["var1"], [1.0])
    @test comp13.eqs_absolute_deviation == Dictionary{EID, Float64}()
  end

  @testset "Compare solutions 4 and 5" begin
    comp45 = compare(sol4, sol5)
    @test comp45.vars_absolute_deviation == Dictionary{Label, Float64}(["var1"], [2.0])
    @test comp45.vars_relative_deviation == Dictionary{Label, Float64}(["var1"], [1.0])
    @test comp45.eqs_absolute_deviation == Dictionary{EID, Float64}([1], [2.0])
  end

  # Check that `compare` between the same model gives zero difference.
  for m in models
    @testset "Compare STAN solution for model $(basename(m)) with itself" begin
      filepath = m
      stan_model = parse(STANTrace, filepath)
      sol_stan = convert(ReconciliationSolution, stan_model)
      stan_variables = keys(sol_stan.variables_value)
      zerodiff_vars = Dictionary{String, Float64}(stan_variables, zeros(Float64, length(stan_variables)))
      stan_equations = keys(sol_stan.score.equations_value)
      zerodiff_eqs = Dictionary{EID, Float64}(stan_equations, zeros(Float64, length(stan_equations)))
      compself = compare(sol_stan, sol_stan)
      @test compself.vars_absolute_deviation == zerodiff_vars
      @test compself.vars_relative_deviation == zerodiff_vars
      @test compself.eqs_absolute_deviation == zerodiff_eqs
    end
  end

  @testset "Compare unobservable variables" begin
    rec = ReconciliationProblem(; name = "Test problem")
    a = MeasuredVariable(; name = "a", value = 10.0, uncertainty = 1.0)
    b = UnmeasuredVariable(; name = "b")
    c = UnmeasuredVariable(; name = "c")
    d = UnmeasuredVariable(; name = "d")
    e = FixedVariable(; name = "e", value = 3.0)
    add_variable!(rec, [a, b, c, d, e])
    rsol1 = ReconciliationSolution(;
      problem = rec,
      variables_value = dictionary(["a" => 11.0, "b" => 5.0, "c" => 1.0, "e" => 3.0]),
      unobservable_variables = dictionary(["c" => c, "d" => d]),
    )
    rsol2 = ReconciliationSolution(;
      problem = rec,
      variables_value = dictionary(["a" => 11.0, "b" => 5.0, "c" => 2.0, "d" => 4.0, "e" => 3.0]),
      unobservable_variables = dictionary(["c" => c, "d" => d]),
    )
    rsol3 = ReconciliationSolution(;
      problem = rec,
      variables_value = dictionary(["a" => 11.0, "b" => 5.0, "c" => 2.0, "d" => 4.0, "e" => 3.0]),
      unobservable_variables = dictionary(["c" => c]),
    )
    rsol4 = ReconciliationSolution(;
      problem = rec,
      variables_value = dictionary(["a" => 11.0, "b" => 5.0, "c" => 2.0, "d" => 4.0, "e" => 3.0]),
      unobservable_variables = dictionary(["b" => b]),
    )

    comp12 = compare(rsol1, rsol2)
    comp23 = compare(rsol2, rsol3)
    comp14 = compare(rsol1, rsol4)

    @test isequal(
      comp12.vars_equal_observability,
      dictionary(["a" => true, "b" => true, "c" => true, "d" => true, "e" => true]),
    )
    @test isequal(
      comp23.vars_equal_observability,
      dictionary(["a" => true, "b" => true, "c" => true, "d" => false, "e" => true]),
    )
    @test isequal(
      comp14.vars_equal_observability,
      dictionary(["a" => true, "b" => false, "c" => false, "d" => false, "e" => true]),
    )
  end
end


@testset "compare jacobian and jacobians_and_hessian in ReconciliationModule" begin
  x = UnmeasuredVariable(; name = "x")
  y = MeasuredVariable(; name = "y", value = 2.0, uncertainty = 1.0)
  z = UnmeasuredVariable(; name = "z")
  v = FixedVariable(; name = "v", value = 3.0)
  lineq1 = Equation(
    ConstantTerm(-6.0),
    LinearTerm[LinearTerm("x", 1.0), LinearTerm("y", 1.0), LinearTerm("z", 1.0)],
    BilinearTerm[],
  )
  lineq2 = Equation(
    ConstantTerm(-8.0),
    LinearTerm[LinearTerm("x", 2.0), LinearTerm("y", 1.0), LinearTerm("z", 1.0)],
    BilinearTerm[],
  )
  lineq3 = Equation(
    ConstantTerm(-14.0),
    LinearTerm[LinearTerm("x", 3.0), LinearTerm("y", 2.0), LinearTerm("z", 2.0)],
    BilinearTerm[],
  ) # lineq3 = lineq1 + lineq2
  bilineq1 = Equation(ConstantTerm(-4.0), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
  bilineq2 = Equation(ConstantTerm(-16.0), LinearTerm[], BilinearTerm[BilinearTerm("x", "z", 2.0)])
  bilineq3 = Equation(
    ConstantTerm(-11.0),
    LinearTerm[LinearTerm("y", 1.0), LinearTerm("v", 1.0)],
    BilinearTerm[BilinearTerm("x", "v", 1.0)],
  )
  problem = ReconciliationProblem(; name = "3 dep vars")
  add_variable!(problem, [x, y, z, v])
  add_equation!(problem, [lineq1, lineq2, lineq3, bilineq1, bilineq2, bilineq3])
  vars_values = Dictionary(["x", "y", "z", "v"], [1.0, 2.0, 3.0, 3.0])
  Jac1 = jacobians_and_hessian(problem, vars_values, zeros(Float64, length(problem.equations)))
  Jac2 = jacobian(problem, vars_values)
  @test Jac1.Jx == Jac2.Jx
  @test Jac1.Jy == Jac2.Jy
end

@testset "compare jacobian and jacobian! in NonlinearSolveInterfaceModule" begin
  x = UnmeasuredVariable(; name = "x")
  y = MeasuredVariable(; name = "y", value = 2.0, uncertainty = 1.0)
  z = UnmeasuredVariable(; name = "z")
  v = FixedVariable(; name = "v", value = 3.0)
  lineq1 = Equation(
    ConstantTerm(-6.0),
    LinearTerm[LinearTerm("x", 1.0), LinearTerm("y", 1.0), LinearTerm("z", 1.0)],
    BilinearTerm[],
  )
  lineq2 = Equation(
    ConstantTerm(-8.0),
    LinearTerm[LinearTerm("x", 2.0), LinearTerm("y", 1.0), LinearTerm("z", 1.0)],
    BilinearTerm[],
  )
  lineq3 = Equation(
    ConstantTerm(-14.0),
    LinearTerm[LinearTerm("x", 3.0), LinearTerm("y", 2.0), LinearTerm("z", 2.0)],
    BilinearTerm[],
  ) # lineq3 = lineq1 + lineq2
  bilineq1 = Equation(ConstantTerm(-4.0), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
  bilineq2 = Equation(ConstantTerm(-16.0), LinearTerm[], BilinearTerm[BilinearTerm("x", "z", 2.0)])
  bilineq3 = Equation(
    ConstantTerm(-11.0),
    LinearTerm[LinearTerm("y", 1.0), LinearTerm("v", 1.0)],
    BilinearTerm[BilinearTerm("x", "v", 1.0)],
  )
  problem = ReconciliationProblem(; name = "3 dep vars")
  add_variable!(problem, [x, y, z, v])
  add_equation!(problem, [lineq1, lineq2, lineq3, bilineq1, bilineq2, bilineq3])
  vars_values = Dictionary(["x", "y", "z", "v"], [1.0, 2.0, 3.0, 3.0])
  problem_qcqp = convert(QCQP, problem)
  J = Array{Float64}(undef, length(problem_qcqp.constraints), length(problem_qcqp.variable_index))
  vars_values_vector = vcat(
    [vars_values[k] for k in keys(problem.measured)],
    [vars_values[k] for k in keys(problem.unmeasured)],
    [vars_values[k] for k in keys(problem.fixed)],
  )
  NonlinearSolveInterfaceModule.jacobian!(J, vars_values_vector, "", problem_qcqp)
  (; Jx, Jy, Jz, b) = jacobian(problem, vars_values)
  @test J == hcat(Jx, Jy, Jz)
end

@testset "Relative residuals" begin
  @testset "compute_relative_residual basic cases" begin
    # Linear equation: 2x + 3y - 5 = 0, at x=2, y=1: residual = 2
    eq_linear = Equation(ConstantTerm(-5.0), LinearTerm[LinearTerm("x", 2.0), LinearTerm("y", 3.0)], BilinearTerm[])
    vars = Dictionary{Label, Float64}(["x", "y"], [2.0, 1.0])
    # denominator = |c| + ||q|| ||x|| = 5 + sqrt(13)*sqrt(5)
    @test ReconciliationModule.compute_relative_residual(eq_linear, 2.0, vars) ≈ 2.0 / (5.0 + sqrt(13.0) * sqrt(5.0)) rtol =
      1e-10

    # Bilinear equation: xy - 6 = 0, at x=3, y=3: residual = 3
    eq_bilinear = Equation(ConstantTerm(-6.0), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    vars2 = Dictionary{Label, Float64}(["x", "y"], [3.0, 3.0])
    # denominator = 6 + 1*18 = 24
    @test ReconciliationModule.compute_relative_residual(eq_bilinear, 3.0, vars2) ≈ 3.0 / 24.0 rtol = 1e-10

    # Zero residual case
    @test ReconciliationModule.compute_relative_residual(eq_linear, 0.0, vars) == 0.0

    # Zero denominator edge case (fallback to absolute residual)
    eq_zero = Equation(ConstantTerm(0.0), LinearTerm[], BilinearTerm[])
    @test ReconciliationModule.compute_relative_residual(eq_zero, 1.0, Dictionary{Label, Float64}()) == 1.0
  end

  @testset "Integration with SolutionScore" begin
    problem = ReconciliationProblem(; name = "Test")
    add_variable!(problem, [UnmeasuredVariable(; name = "x"), UnmeasuredVariable(; name = "y")])
    add_equation!(
      problem,
      Equation(ConstantTerm(-10.0), LinearTerm[LinearTerm("x", 3.0), LinearTerm("y", 4.0)], BilinearTerm[]),
    )

    vars_value = Dictionary{Label, Float64}(["x", "y"], [1.0, 1.0])
    score = ReconciliationModule.evaluate_score(problem, vars_value)

    # Check relative_residuals field exists and is computed
    @test haskey(score.relative_residuals, 1)
    @test abs(score.equations_value[1]) == 3.0
    @test score.relative_residuals[1] ≈ 3.0 / (10.0 + 5.0 * sqrt(2.0)) rtol = 1e-10
  end
end

@testset "Computation of the Jacobian" begin
  include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))
  problem = seven_flows_cencic_2012()
  remove_fixed_variables!(problem)
  vals = dictionary([
    # Measured variables
    "m1" => 102.42595948021597,
    "m3" => 302.41622972562095,
    "m5" => 152.42595998014826,
    "tc34" => 0.4959729505319356,
    # Unmeasured variables
    "m4" => 149.99026974571603,
    "m6" => 76.21297999007412,
    "m7" => 76.21297999007415,
  ])
  @test norm(evaluate(problem, vals), Inf) <= 1e-6

  (; Jx, Jy, variable_index) = jacobians(problem, vals)
  #= Problem equations:
   1 │  m1 + m4 - m3 = -50.0
   2 │  m3 - m4 - m5 = 0
   3 │  m5 - m6 - m7 = 0
   4 │  m4 - m3 * tc34 = 0
  =#
  @test issparse(Jx)
  @test Jx == [1 -1 0 0; 0 1 -1 0; 0 0 1 0; 0 -vals["tc34"] 0 -vals["m3"]]
  @test issparse(Jy)
  @test Jy == [1 0 0; -1 0 0; 0 -1 -1; 1 0 0]
  # This problem doesn't have full column rank. The solution space at the linearization
  # point has dimension 1.
  @test rank(Jy) == 2
  @test rref(Jy) == [1 0 0; 0 1 1; 0 0 0; 0 0 0]
end

@testset "Status evaluation" begin
  if length(all_models()) >= 17
    m = all_models()[17]
    strace = parse(STANTrace, m)
    sol_stan = convert(ReconciliationSolution, strace)
    @test isequal(sol_stan.status, FEASIBLE_WITH_UNOBSERVABLE_VARIABLES)
  else
    @test_skip "models directory not available (expected >= 17 models)"
  end
end
