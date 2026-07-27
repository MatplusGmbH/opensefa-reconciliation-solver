using Test
using ResultTypes: iserror, unwrap
using Dictionaries
using LinearAlgebra: rank
using SparseArrays
import OpenSEFA
using OpenSEFA.STANModule
using OpenSEFA.ConstraintsModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.PresolveModule
using OpenSEFA.ErrorModule
using OpenSEFA.ErrorInfosModule
using OpenSEFA.ConstraintParserModule
using OpenSEFA.DefaultSolverModule
using OpenSEFA.ComparisonModule
using OpenSEFA.SerializationModule
using OpenSEFA.JuMPInterfaceModule

# Filepath for each STAN model available under the `models` folder.
models = all_models()

@testset "replace variables" begin
  @testset "replace variables including square term" begin
    x = UnmeasuredVariable(; name = "x")
    y = UnmeasuredVariable(; name = "y")
    eq = Equation(ConstantTerm(0.0), [LinearTerm("y", 1.0)], [BilinearTerm("x", "x", 1.0)])
    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)
    replace_variables!(problem, Dictionary(["x"], [2.0]))
    @test problem.equations[1] == Equation(ConstantTerm(4.0), [LinearTerm("y", 1.0)], BilinearTerm[])
  end

  @testset "replace variables such that linear terms cancel" begin
    x = UnmeasuredVariable(; name = "x")
    y = UnmeasuredVariable(; name = "y")
    eq = Equation(ConstantTerm(3.0), [LinearTerm("x", -4.0), LinearTerm("y", 2.0)], [BilinearTerm("x", "y", 2.0)])
    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)
    replace_variables!(problem, Dictionary(["y"], [2.0]))
    @test problem.equations[1] == Equation(ConstantTerm(7.0), LinearTerm[], BilinearTerm[])
  end

  @testset "replace variables resulting in variable removal" begin
    x = UnmeasuredVariable(; name = "x")
    y = UnmeasuredVariable(; name = "y")
    z = UnmeasuredVariable(; name = "z")
    v = UnmeasuredVariable(; name = "v")
    w = UnmeasuredVariable(; name = "w")
    eq = Equation(ConstantTerm(0.0), [LinearTerm("x", -4.0)], [BilinearTerm("x", "y", 2.0)])
    eq2 = Equation(ConstantTerm(2.0), [LinearTerm("x", 2.0)], [BilinearTerm("v", "w", 3.0)])
    eq3 = Equation(ConstantTerm(3.0), [LinearTerm("x", -4.0), LinearTerm("z", 1.0)], [BilinearTerm("x", "y", 2.0)])
    problem = ReconciliationProblem()
    add_variable!(problem, [x, y, z, v, w])
    add_equation!(problem, [eq, eq2, eq3])
    presolved = unwrap(replace_variables!(problem, Dictionary(["y", "v"], [2.0, 0.0])))
    @test collect(keys(presolved.pre.variables)) == ["x", "z"]
    @test !(1 in keys(presolved.pre.equations))
    @test presolved.pre.equations[2] == Equation(ConstantTerm(2.0), [LinearTerm("x", 2.0)], BilinearTerm[])
    @test presolved.pre.equations[3] == Equation(ConstantTerm(3.0), [LinearTerm("z", 1.0)], BilinearTerm[])
    @test presolved.pre.var_to_eq == Dictionary(["x", "z"], [[2], [3]])
    @test presolved.pre.eq_to_var == Dictionary([2, 3], [["x"], ["z"]])
    @test collect(keys(presolved.unobservable)) == ["w"]
  end
end

@testset "Replacement of products." begin
  @testset "One equation without storing equation" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    eq = Equation(
      LinearTerm[LinearTerm("x", 1.0)],
      BilinearTerm[BilinearTerm("x", "y", 1.0)])
    rec = ReconciliationProblem(name = "One equation with product")
    add_variable!(rec, x)
    add_variable!(rec, y)
    add_equation!(rec, eq)
    replace_bilinear_terms!(rec, Dictionary{Tuple{Label, Label}, Float64}([(x.name, y.name)], [1.0]))
    @test rec.equations[1].linear_terms == LinearTerm[LinearTerm("x", 1.0)]
    @test rec.equations[1].bilinear_terms == []
    @test rec.equations[1].constant_term == ConstantTerm(1)
    @test rec.var_to_eq == Dictionary{Label, Vector{EID}}([x.name, y.name], [[1], []])
    @test rec.eq_to_var == Dictionary{EID, Vector{Label}}([1], [[x.name]])
  end
  @testset "Two equations with permutation without storing equation" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    eq1 = Equation(
      LinearTerm[LinearTerm("x", 1.0)],
      BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(
      LinearTerm[LinearTerm("x", 2.0)],
      BilinearTerm[BilinearTerm("y", "x", 2.0)])
    rec = ReconciliationProblem(name = "Two equations with permutated products")
    add_variable!(rec, x)
    add_variable!(rec, y)
    add_equation!(rec, [eq1, eq2])
    replace_bilinear_terms!(
      rec,
      Dictionary{Tuple{Label, Label}, Float64}([(x.name, y.name), (y.name, x.name)], [1.0, 1.0]),
      Dictionary{Tuple{Label, Label}, EID}(),
    )
    @test rec.equations[1].linear_terms == LinearTerm[LinearTerm("x", 1.0)]
    @test rec.equations[1].bilinear_terms == []
    @test rec.equations[1].constant_term == ConstantTerm(1)
    @test rec.equations[2].linear_terms == LinearTerm[LinearTerm("x", 2.0)]
    @test rec.equations[2].bilinear_terms == []
    @test rec.equations[2].constant_term == ConstantTerm(2)
    @test rec.var_to_eq == Dictionary{Label, Vector{EID}}([x.name, y.name], [[1, 2], []])
    @test rec.eq_to_var == Dictionary{EID, Vector{Label}}([1, 2], [[x.name], [x.name]])
  end
  @testset "Two equations with several bilinear terms without storing equation" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(
      LinearTerm[LinearTerm("x", 1.0)],
      BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(
      LinearTerm[LinearTerm("x", 2.0)],
      BilinearTerm[BilinearTerm("x", "y", 1.0), BilinearTerm("y", "z", 2.0)])
    rec = ReconciliationProblem(name = "Two equations several bilinear terms")
    add_variable!(rec, [x, y, z])
    add_equation!(rec, [eq1, eq2])
    replace_bilinear_terms!(
      rec,
      Dictionary{Tuple{Label, Label}, Float64}([(x.name, y.name)], [1.0]),
      Dictionary{Tuple{Label, Label}, EID}(),
    )
    @test rec.equations[1].linear_terms == LinearTerm[LinearTerm("x", 1.0)]
    @test rec.equations[1].bilinear_terms == []
    @test rec.equations[1].constant_term == ConstantTerm(1)
    @test rec.equations[2].linear_terms == LinearTerm[LinearTerm("x", 2.0)]
    @test rec.equations[2].bilinear_terms == [BilinearTerm("y", "z", 2.0)]
    @test rec.equations[2].constant_term == ConstantTerm(1)
    @test rec.var_to_eq == Dictionary{Label, Vector{EID}}([x.name, y.name, z.name], [[1, 2], [2], [2]])
    @test rec.eq_to_var == Dictionary{EID, Vector{Label}}([1, 2], [[x.name], [x.name, y.name, z.name]])
  end
  @testset "Three equations with storing equation" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(
      LinearTerm[LinearTerm("x", 1.0)],
      BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-2), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 2.0)])
    eq2_backup = deepcopy(eq2)
    eq3 = Equation(
      LinearTerm[LinearTerm("y", 1.0)],
      BilinearTerm[BilinearTerm("x", "y", 1.0), BilinearTerm("y", "z", 1.0)])
    rec = ReconciliationProblem(name = "Three equations with several bilinear terms")
    add_variable!(rec, [x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    replace_bilinear_terms!(
      rec,
      Dictionary{Tuple{Label, Label}, Float64}([(x.name, y.name), (y.name, x.name)], [1.0, 1.0]),
      Dictionary{Tuple{Label, Label}, EID}([(x.name, y.name)], [2]),
    )
    @test rec.equations[1].linear_terms == LinearTerm[LinearTerm("x", 1.0)]
    @test rec.equations[1].bilinear_terms == []
    @test rec.equations[1].constant_term == ConstantTerm(1)
    @test eq2 == eq2_backup
    @test rec.equations[3].linear_terms == LinearTerm[LinearTerm("y", 1.0)]
    @test rec.equations[3].bilinear_terms == [BilinearTerm("y", "z", 1.0)]
    @test rec.equations[3].constant_term == ConstantTerm(1)
    @test rec.var_to_eq == Dictionary{Label, Vector{EID}}([x.name, y.name, z.name], [[1, 2], [2, 3], [3]])
    @test rec.eq_to_var == Dictionary{EID, Vector{Label}}([1, 2, 3], [[x.name], [x.name, y.name], [y.name, z.name]])
  end
  @testset "Conflicting replacement" begin
    x = UnmeasuredVariable(name = "x")
    y = UnmeasuredVariable(name = "y")
    eq1 = Equation(ConstantTerm(2), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    rec = ReconciliationProblem(name = "Three equations with several bilinear terms")
    add_variable!(rec, [x, y])
    add_equation!(rec, [eq1])
    res = replace_bilinear_terms!(
      rec,
      Dictionary{Tuple{Label, Label}, Float64}([(x.name, y.name), (y.name, x.name)], [1.0, 1.0]),
    )
    @test iserror(res)
    @test_throws InfeasibleEquationError unwrap(res)
  end
end

@testset "Removal of fixed variables." begin
  @testset "Single bilinear term." begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    z1 = FixedVariable(name = "z1", value = 8.0)
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0)],
      BilinearTerm[BilinearTerm("z1", "z1", 2.0)])
    rec = ReconciliationProblem(name = "Single bilinear term.")
    add_variable!(rec, x1)
    add_variable!(rec, z1)
    add_equation!(rec, eq1)
    remove_fixed_variables!(rec)
    @test isempty(rec.fixed)
    @test rec.equations[1] ==
          Equation(ConstantTerm(128.0), LinearTerm[LinearTerm("x1", 1.0)], BilinearTerm[])
  end

  @testset "One-dimensional case." begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5)
    z1 = FixedVariable(name = "z1", value = 8.0)
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("z1", -1.0)],
      BilinearTerm[])
    eq2 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x2", 1.0)],
      BilinearTerm[BilinearTerm("x1", "z1", 1.0)])
    rec = ReconciliationProblem(name = "One fixed variable")
    add_variable!(rec, x1)
    add_variable!(rec, x2)
    add_variable!(rec, z1)
    add_equation!(rec, [eq1, eq2])
    remove_fixed_variables!(rec)
    @test isempty(rec.fixed)
    @test rec.equations[1] ==
          Equation(ConstantTerm(-8.0), LinearTerm[LinearTerm("x1", 1.0)], BilinearTerm[])
    @test rec.equations[2] ==
          Equation(LinearTerm[LinearTerm("x1", 9.0), LinearTerm("x2", 1.0)], BilinearTerm[])
  end

  @testset "Integration tests." begin
    for m in models
      # Convert data to reconciliation problem.
      trace_filename = joinpath(m, "data", "trace.txt")
      if !isfile(trace_filename)
        @warn "Skipping model $(basename(m)). Expected trace file under `/data/trace.txt`."
        continue
      end
      strace = parse(STANTrace, m)
      rec = convert(ReconciliationProblem, strace)

      remove_fixed_variables!(rec)
      @test isempty(rec.fixed)
    end
  end
end

@testset "Removal of redundant equations" begin
  @testset "Linearly dependent equations constant case" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)

    # Equation in the form of -1 + x1 = 0
    eq1 = Equation(
      ConstantTerm(-1.0),
      LinearTerm[LinearTerm("x1", 1.0)],
      BilinearTerm[])

    # Equation in the form of -2 + 2x1 = 0
    eq2 = Equation(
      ConstantTerm(-2.0),
      LinearTerm[LinearTerm("x1", 2.0)],
      BilinearTerm[])

    rec = ReconciliationProblem(name = "Linear with Constant Case")
    add_variable!(rec, x1)
    add_equation!(rec, [eq1, eq2])
    remove_redundant_equations!(rec)

    # Test for the removal of eq2 and the remaining equation eq1
    @test only(rec.equations) == eq1
  end

  @testset "Linearly dependent equations linear case" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5)

    # Equation in the form of x1 + x2 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x2", 1.0)],
      BilinearTerm[])

    # Equation in the form of 2x1 + 2x2 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("x1", 2.0), LinearTerm("x2", 2.0)],
      BilinearTerm[])

    rec = ReconciliationProblem(name = "Linear Case")
    add_variable!(rec, [x1, x2])
    add_equation!(rec, [eq1, eq2])
    remove_redundant_equations!(rec)

    # Test for the removal of eq2 and the remaining equation eq1
    @test only(rec.equations) == eq1
  end

  @testset "Linearly dependent equations bilinear case" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5)

    # Equation in the form of x1x2 + x1 + x2 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x2", 1.0)],
      BilinearTerm[BilinearTerm("x1", "x2", 1.0)])

    # Equation in the form of 2x1x2 + 2x1 + 2x2 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("x1", 2.0), LinearTerm("x2", 2.0)],
      BilinearTerm[BilinearTerm("x1", "x2", 2.0)])

    rec = ReconciliationProblem(name = "Bilinear and Linear Case")
    add_variable!(rec, [x1, x2])
    add_equation!(rec, [eq1, eq2])
    remove_redundant_equations!(rec)

    # Test for the removal of eq2 and the remaining equation eq1
    @test only(rec.equations) == eq1
  end
end

@testset "Canonical Replacement" begin
  @testset "linear equations" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5, lb = -Inf)
    z1 = FixedVariable(name = "z1", value = 8.0)

    # Equation in the form of x1 - z1 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("z1", -1.0)],
      BilinearTerm[])

    # Equation in the form of x1 + x2 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x2", 1.0)],
      BilinearTerm[])

    rec = ReconciliationProblem(name = "Linear replacement")
    add_variable!(rec, [z1, x1, x2])
    add_equation!(rec, [eq1, eq2])
    presolved = unwrap(remove_fixed_variables!(rec))
    result = trivial_substitutions_recursive!(presolved)

    @test !iserror(result)
    @test length(presolved.pre.equations) == 0
    @test length(presolved.replacements) == 3

    # Type of replacements is Dictionaries.Dictionary not Dict of Julia so it is converted into Dict
    @test Dict(pairs(presolved.replacements)) == Dict("z1" => 8.0, "x1" => 8.0, "x2" => -8.0)
  end

  @testset "bilinear equations" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5)
    z1 = FixedVariable(name = "z1", value = 2.0)

    # Equation in the form of x1 - z1 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("z1", -1.0)],
      BilinearTerm[])

    # Equation in the form of x1x2 + x2 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("x2", 1.0)],
      BilinearTerm[BilinearTerm("x1", "x2", 1.0)])

    rec = ReconciliationProblem(name = "Bilinear replacement")
    add_variable!(rec, [z1, x1, x2])
    add_equation!(rec, [eq1, eq2])

    presolved = unwrap(remove_fixed_variables!(rec))
    result = trivial_substitutions_recursive!(presolved)

    @test !iserror(result)
    @test length(presolved.pre.equations) == 0
    @test length(presolved.replacements) == 3
    @test Dict(pairs(presolved.replacements)) == Dict("z1" => 2.0, "x1" => 2.0, "x2" => 0.0)
  end

  @testset "bilinear and linear terms" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5, lb = -Inf)
    z1 = FixedVariable(name = "z1", value = 2.0)

    # Equation in the form of x1 - z1 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("z1", -1.0)],
      BilinearTerm[])

    # Equation in the form of x1x2 + x1 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("x1", 1.0)],
      BilinearTerm[BilinearTerm("x1", "x2", 1.0)])

    rec = ReconciliationProblem(name = "Bilinear and linear replacement")
    add_variable!(rec, [z1, x1, x2])
    add_equation!(rec, [eq1, eq2])

    presolved = unwrap(remove_fixed_variables!(rec))
    result = trivial_substitutions_recursive!(presolved)

    @test !iserror(result)
    @test length(presolved.pre.equations) == 0
    @test length(presolved.replacements) == 3
    @test Dict(pairs(presolved.replacements)) == Dict("z1" => 2.0, "x1" => 2.0, "x2" => -1.0)
  end

  @testset "square terms" begin
    x = FixedVariable(; name = "x", value = 2.0)
    y = UnmeasuredVariable(; name = "y", lb = -Inf)
    eq = parse_equation("x * x + y = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)

    pre = presolve(problem)
    replacements = unwrap(pre).replacements

    @test Dict(pairs(replacements)) == Dict("x" => 2.0, "y" => -4.0)
  end

  @testset "square replacements" begin
    x = FixedVariable(; name = "x", value = -4.0)
    y = UnmeasuredVariable(; name = "y")
    eq = parse_equation("x + y * y = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)

    pre = presolve(problem)
    replacements = unwrap(pre).replacements

    @test Dict(pairs(replacements)) == Dict("x" => -4.0, "y" => 2.0)
  end

  @testset "square replacements (negative square root solution)" begin
    x = FixedVariable(; name = "x", value = -4.0)
    y = UnmeasuredVariable(; name = "y", lb = -10.0, ub = 1.0)
    eq = parse_equation("x + y * y = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)

    pre = presolve(problem)
    replacements = unwrap(pre).replacements

    @test Dict(pairs(replacements)) == Dict("x" => -4.0, "y" => -2.0)
  end

  @testset "square replacements (no unique square root solution)" begin
    x = FixedVariable(; name = "x", value = -4.0)
    y = UnmeasuredVariable(; name = "y", lb = -10.0)
    eq = parse_equation("x + y * y = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)

    pre = presolve(problem)
    replacements = unwrap(pre).replacements

    @test Dict(pairs(replacements)) == Dict("x" => -4.0)
  end

  @testset "square replacements (value 0.0)" begin
    y = UnmeasuredVariable(; name = "y")
    eq = parse_equation("y * y = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [y])
    add_equation!(problem, eq)

    pre = presolve(problem)
    replacements = unwrap(pre).replacements

    @test Dict(pairs(replacements)) == Dict("y" => 0.0)
  end

  @testset "square replacements (invalid: negative square)" begin
    x = FixedVariable(; name = "x", value = 4.0)
    y = UnmeasuredVariable(; name = "y")
    eq = parse_equation("x + y * y = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)

    pre = presolve(problem)

    @test iserror(pre)
  end

  @testset "square replacements (invalid: no solution within variable bounds)" begin
    x = FixedVariable(; name = "x", value = -4.0)
    y = UnmeasuredVariable(; name = "y", ub = 1.0)
    eq = parse_equation("x + y * y = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [x, y])
    add_equation!(problem, eq)

    pre = presolve(problem)

    @test iserror(pre)
  end
end

@testset "Robustness to small factors in substitutions" begin
  # Regression test for issue #108: a near-zero factor such as 1e-15 should not trigger
  # a substitution. The equation `1e-15 * x = 0` is essentially `0 = 0` due to rounding.
  @testset "near-zero factor" begin
    x = UnmeasuredVariable(; name = "x")
    eq = parse_equation("0.000000000000001 * x = 0")
    problem = ReconciliationProblem()
    add_variable!(problem, x)
    add_equation!(problem, eq)

    prerec = PresolvedReconciliationProblem(problem)
    result = trivial_substitutions!(prerec)
    @test !iserror(result)
    # No substitution should be performed: x remains free
    @test isempty(prerec.replacements)
  end

  # A non-near-zero factor with a near-zero constant should substitute x = 0,
  # not x = -c/t ≈ some spurious non-zero value from floating-point rounding.
  @testset "near-zero constant term normalised before division" begin
    x = UnmeasuredVariable(; name = "x")
    # 2.0 * x = 1e-15: the RHS is near-zero, so x should be substituted as 0.
    rhs = 1e-15  # isapproxzero(rhs) is true
    eq = Equation(ConstantTerm(-rhs), LinearTerm[LinearTerm("x", 2.0)], BilinearTerm[])
    problem = ReconciliationProblem()
    add_variable!(problem, x)
    add_equation!(problem, eq)

    prerec = PresolvedReconciliationProblem(problem)
    result = trivial_substitutions!(prerec)
    @test !iserror(result)
    @test Dict(pairs(prerec.replacements)) == Dict("x" => 0.0)
  end
end

@testset "Bilinear product replacements" begin
  @testset "Three equations" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(LinearTerm[LinearTerm("x", 1.0)], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-2), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 2.0)])
    eq2_backup = deepcopy(eq2)
    eq3 =
      Equation(LinearTerm[LinearTerm("y", 1.0)], BilinearTerm[BilinearTerm("x", "y", 1.0), BilinearTerm("y", "z", 1.0)])
    rec = ReconciliationProblem(name = "Three equations")
    add_variable!(rec, [x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    bilinear_substitutions!(rec)
    @test eq1 == Equation(ConstantTerm(1), LinearTerm[LinearTerm("x", 1.0)], BilinearTerm[])
    @test eq2 == eq2_backup
    @test eq3 == Equation(ConstantTerm(1), LinearTerm[LinearTerm("y", 1.0)], BilinearTerm[BilinearTerm("y", "z", 1.0)])
  end
  @testset "Three equations with permutation" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(LinearTerm[LinearTerm("x", 1.0)], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-2), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 2.0)])
    eq2_backup = deepcopy(eq2)
    eq3 =
      Equation(LinearTerm[LinearTerm("y", 1.0)], BilinearTerm[BilinearTerm("y", "x", 1.0), BilinearTerm("y", "z", 1.0)])
    rec = ReconciliationProblem(name = "Three equations with permutation")
    add_variable!(rec, [x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    bilinear_substitutions!(rec)
    @test eq1 == Equation(ConstantTerm(1), LinearTerm[LinearTerm("x", 1.0)], BilinearTerm[])
    @test eq2 == eq2_backup
    @test eq3 == Equation(ConstantTerm(1), LinearTerm[LinearTerm("y", 1.0)], BilinearTerm[BilinearTerm("y", "z", 1.0)])
  end
  @testset "Two replacements in one iteration" begin
    v = UnmeasuredVariable(name = "v")
    w = UnmeasuredVariable(name = "w")
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(ConstantTerm(-1), LinearTerm[], BilinearTerm[BilinearTerm("v", "w", 1.0)])
    eq1_backup = deepcopy(eq1)
    eq2 = Equation(ConstantTerm(-6), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 2.0)])
    eq2_backup = deepcopy(eq2)
    eq3 = Equation(
      ConstantTerm(-7),
      LinearTerm[LinearTerm("z", 1.0)],
      BilinearTerm[BilinearTerm("v", "w", 1.0), BilinearTerm("y", "x", 1.0)],
    )
    rec = ReconciliationProblem(name = "Three equations")
    add_variable!(rec, [v, w, x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    bilinear_substitutions!(rec)
    @test eq1 == eq1_backup
    @test eq2 == eq2_backup
    @test eq3 == Equation(ConstantTerm(-3), LinearTerm[LinearTerm("z", 1.0)], BilinearTerm[])
  end
  @testset "Two replacements in subsequent iterations" begin
    v = UnmeasuredVariable(name = "v")
    w = UnmeasuredVariable(name = "w")
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 =
      Equation(ConstantTerm(-6), LinearTerm[], BilinearTerm[BilinearTerm("w", "v", 2.0), BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-1), LinearTerm[], BilinearTerm[BilinearTerm("v", "w", 1.0)])
    eq2_backup = deepcopy(eq2)
    eq3 = Equation(
      ConstantTerm(-7),
      LinearTerm[LinearTerm("z", 1.0)],
      BilinearTerm[BilinearTerm("v", "w", 1.0), BilinearTerm("y", "x", 1.0)],
    )
    rec = ReconciliationProblem(name = "Three equations")
    add_variable!(rec, [v, w, x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    bilinear_substitutions!(rec)
    @test eq1 == Equation(ConstantTerm(-4), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    @test eq2 == eq2_backup
    @test eq3 == Equation(ConstantTerm(-6), LinearTerm[LinearTerm("z", 1.0)], BilinearTerm[BilinearTerm("y", "x", 1.0)])
    bilinear_substitutions!(rec)
    @test eq1 == Equation(ConstantTerm(-4), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    @test eq2 == eq2_backup
    @test eq3 == Equation(ConstantTerm(-2), LinearTerm[LinearTerm("z", 1.0)], BilinearTerm[])
  end
  @testset "Conflicting equations" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    eq1 = Equation(ConstantTerm(-1), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-6), LinearTerm[], BilinearTerm[BilinearTerm("y", "x", 2.0)])
    rec = ReconciliationProblem(name = "Two conflicting equations")
    add_variable!(rec, [x, y])
    add_equation!(rec, [eq1, eq2])
    res = bilinear_substitutions!(rec)
    @test iserror(res)
    @test_throws InfeasibleEquationError unwrap(res)
  end
end

@testset "Recursive trivial and product replacements" begin
  @testset "Three equations" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(LinearTerm[LinearTerm("x", 1.0)], BilinearTerm[BilinearTerm("x", "y", -1.0)])
    eq2 = Equation(ConstantTerm(-2), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 2.0)])
    eq2_backup = deepcopy(eq2)
    eq3 =
      Equation(
        LinearTerm[LinearTerm("y", 1.0)],
        BilinearTerm[BilinearTerm("x", "y", 1.0), BilinearTerm("y", "z", -1.0)],
      )
    rec = ReconciliationProblem(name = "Three equations")
    add_variable!(rec, [x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    res = trivial_and_bilinear_substitutions_recursive!(rec)
    repl = unwrap(res).replacements
    @test isequal(rec, ReconciliationProblem(name = rec.name))
    @test isequal(Dict(pairs(repl)), Dict("x" => 1.0, "y" => 1.0, "z" => 2.0))
  end
  @testset "Three equations with permutation" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(LinearTerm[LinearTerm("x", 1.0)], BilinearTerm[BilinearTerm("x", "y", -1.0)])
    eq2 = Equation(ConstantTerm(-2), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 2.0)])
    eq2_backup = deepcopy(eq2)
    eq3 =
      Equation(
        LinearTerm[LinearTerm("y", 1.0)],
        BilinearTerm[BilinearTerm("y", "x", 1.0), BilinearTerm("y", "z", -1.0)],
      )
    rec = ReconciliationProblem(name = "Three equations with permutation")
    add_variable!(rec, [x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    res = trivial_and_bilinear_substitutions_recursive!(rec)
    repl = unwrap(res).replacements
    @test isequal(rec, ReconciliationProblem(name = rec.name))
    @test isequal(Dict(pairs(repl)), Dict("x" => 1.0, "y" => 1.0, "z" => 2.0))
  end
  @testset "Two replacements in one iteration" begin
    v = UnmeasuredVariable(name = "v")
    w = UnmeasuredVariable(name = "w")
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(ConstantTerm(-1), LinearTerm[], BilinearTerm[BilinearTerm("v", "w", 1.0)])
    eq1_backup = deepcopy(eq1)
    eq2 = Equation(ConstantTerm(-6), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 2.0)])
    eq2_backup = deepcopy(eq2)
    eq3 = Equation(
      ConstantTerm(-7),
      LinearTerm[LinearTerm("z", 1.0)],
      BilinearTerm[BilinearTerm("v", "w", 1.0), BilinearTerm("y", "x", 1.0)],
    )
    rec = ReconciliationProblem(name = "Three equations")
    add_variable!(rec, [v, w, x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    res = trivial_and_bilinear_substitutions_recursive!(rec)
    repl = unwrap(res).replacements
    @test eq1 == eq1_backup
    @test eq2 == eq2_backup
    @test !haskey(rec.equations, 3)
    @test !haskey(rec.variables, z.name)
    @test isequal(Dict(pairs(repl)), Dict("z" => 3.0))
  end
  @testset "Two replacements in subsequent iterations" begin
    v = UnmeasuredVariable(name = "v")
    w = UnmeasuredVariable(name = "w")
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 =
      Equation(ConstantTerm(-6), LinearTerm[], BilinearTerm[BilinearTerm("w", "v", 2.0), BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-1), LinearTerm[], BilinearTerm[BilinearTerm("v", "w", 1.0)])
    eq2_backup = deepcopy(eq2)
    eq3 = Equation(
      ConstantTerm(-7),
      LinearTerm[LinearTerm("z", 1.0)],
      BilinearTerm[BilinearTerm("v", "w", 1.0), BilinearTerm("y", "x", 1.0)],
    )
    rec = ReconciliationProblem(name = "Three equations")
    add_variable!(rec, [v, w, x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    res = trivial_and_bilinear_substitutions_recursive!(rec)
    repl = unwrap(res).replacements
    @test eq1 == Equation(ConstantTerm(-4), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    @test eq2 == eq2_backup
    @test !haskey(rec.equations, 3)
    @test !haskey(rec.variables, z.name)
    @test isequal(Dict(pairs(repl)), Dict("z" => 2.0))
  end
  @testset "Conflicting bilinear equations" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    eq1 = Equation(ConstantTerm(-1), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-6), LinearTerm[], BilinearTerm[BilinearTerm("y", "x", 2.0)])
    rec = ReconciliationProblem(name = "Two conflicting equations")
    add_variable!(rec, [x, y])
    add_equation!(rec, [eq1, eq2])
    res = trivial_and_bilinear_substitutions_recursive!(rec)
    @test iserror(res)
    @test_throws InfeasibleEquationError unwrap(res)
    @test isequal(
      Dict(pairs(res.error.infos[INFO_CONFLICTING_REPLACEMENTS])),
      Dict(("x", "y") => [1.0, 3.0], ("y", "x") => [1.0, 3.0]),
    )
    @test isempty(res.error.infos[INFO_PREVIOUS_REPLACEMENTS])
    @test isempty(res.error.infos["Previous product replacements"])
  end
  @testset "Conflicting linear and bilinear equation" begin
    x = MeasuredVariable(name = "x", value = 10.0, uncertainty = 1.0)
    y = UnmeasuredVariable(name = "y")
    z = UnmeasuredVariable(name = "z")
    eq1 = Equation(ConstantTerm(-1), LinearTerm[], BilinearTerm[BilinearTerm("x", "y", 1.0)])
    eq2 = Equation(ConstantTerm(-6), LinearTerm[LinearTerm("z", 2.0)], BilinearTerm[BilinearTerm("y", "x", 2.0)])
    eq3 = Equation(ConstantTerm(-6), LinearTerm[LinearTerm("z", 1.0)], BilinearTerm[BilinearTerm("x", "y", 3.0)])
    rec = ReconciliationProblem(name = "Three equations with conflict")
    add_variable!(rec, [x, y, z])
    add_equation!(rec, [eq1, eq2, eq3])
    res = trivial_and_bilinear_substitutions_recursive!(rec)
    @test iserror(res)
    @test_throws InfeasibleEquationError unwrap(res)
    @test isequal(Dict(pairs(res.error.infos[INFO_CONFLICTING_REPLACEMENTS])), Dict("z" => [2.0, 3.0]))
    @test isempty(res.error.infos[INFO_PREVIOUS_REPLACEMENTS])
    @test isequal(
      Dict(pairs(res.error.infos["Previous product replacements"])),
      Dict(("x", "y") => 1.0, ("y", "x") => 1.0),
    )
  end
end

@testset "Infeasible Cases" begin
  @testset "Spurious equation" begin
    # Let z1 - z2 = 0 => z1 = z2, but z1 = 2.0 and z2 = 3.0, hence this equation is infeasible.
    # Moreover, upon substitution the residual should be -1.0
    z1 = FixedVariable(name = "z1", value = 2.0)
    z2 = FixedVariable(name = "z2", value = 3.0)
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 1.0)

    # Equation in the form of z1 - z2 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("z1", 1.0), LinearTerm("z2", -1.0)],
      BilinearTerm[],
    )

    # Equation in the form of z2 - x1 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("z2", 1.0), LinearTerm("x1", -1.0)],
      BilinearTerm[],
    )

    # Equation in the form of z1 - x2 = 0
    eq3 = Equation(
      LinearTerm[LinearTerm("z1", 1.0), LinearTerm("x2", -1.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Spurious equation")
    add_variable!(rec, z1)
    add_variable!(rec, z2)
    add_variable!(rec, x1)
    add_variable!(rec, x2)
    add_equation!(rec, [eq1, eq2, eq3])

    res = presolve(rec)
    @test iserror(res)
    @test res.error isa InfeasibleEquationError
    @test res.error.index == [1]
    @test res.error.residual == [-1.0]
  end

  @testset "Inconsistent equations" begin
    z1 = FixedVariable(name = "z1", value = 2.0)
    z2 = FixedVariable(name = "z2", value = 3.0)
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5)

    # Equation in the form: - z1 + x1 + x2 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("z1", -1.0), LinearTerm("x1", 1.0), LinearTerm("x2", 1.0)],
      BilinearTerm[],
    )

    # Equation in the form: - z2 + x1 + x2 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("z2", -1.0), LinearTerm("x1", 1.0), LinearTerm("x2", 1.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Inconsistent equations")
    add_variable!(rec, z1)
    add_variable!(rec, x1)
    add_variable!(rec, x2)
    add_variable!(rec, z2)
    add_equation!(rec, [eq1, eq2])

    res = presolve(rec)
    @test iserror(res)
    @test issetequal(keys(res.error.infos[INFO_CONFLICTING_EQ_DEPS]), [2])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][2], [1, 2])
  end

  @testset "Circular Dependencies" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = MeasuredVariable(name = "x2", value = 5.0, uncertainty = 0.5)
    x3 = MeasuredVariable(name = "x3", value = 3.0, uncertainty = 0.3)

    # x1 - x2 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x2", -1.0)],
      BilinearTerm[],
    )

    # x2 - x3 = 0
    eq2 = Equation(
      LinearTerm[LinearTerm("x2", 1.0), LinearTerm("x3", -1.0)],
      BilinearTerm[],
    )

    # x3 + x1 = 0
    eq3 = Equation(
      LinearTerm[LinearTerm("x3", 1.0), LinearTerm("x1", 1.0)],
      BilinearTerm[],
    )

    # x1 + x2 + x3 = 1
    eq4 = Equation(
      ConstantTerm(1.0),
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x2", 1.0), LinearTerm("x3", 1.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Circular Dependencies")
    add_variable!(rec, [x1, x2, x3])
    add_equation!(rec, [eq1, eq2, eq3, eq4])

    res = presolve(rec)
    @test iserror(res)
    @test issetequal(keys(res.error.infos[INFO_CONFLICTING_EQ_DEPS]), [4])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][4], [1, 2, 3, 4])
  end

  @testset "Mixing consistent and inconsistent equations" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = UnmeasuredVariable(name = "x2")
    x3 = MeasuredVariable(name = "x3", value = 3.0, uncertainty = 0.3)
    x4 = MeasuredVariable(name = "x4", value = 5.0, uncertainty = 0.5)

    # x1 - x4 = 0
    eq1 = Equation(
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x4", -1.0)],
      BilinearTerm[],
    )

    # x3 = 3
    eq2 = Equation(
      ConstantTerm(-3.0),
      LinearTerm[LinearTerm("x3", 1.0)],
      BilinearTerm[],
    )

    # x1 + x4 = 2
    eq3 = Equation(
      ConstantTerm(-2.0),
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x4", 1.0)],
      BilinearTerm[],
    )

    # 2*x1 + 3*x4 = 4
    eq4 = Equation(
      ConstantTerm(-4.0),
      LinearTerm[LinearTerm("x1", 2.0), LinearTerm("x4", 3.0)],
      BilinearTerm[],
    )

    # x2 = 2
    eq5 = Equation(
      ConstantTerm(-2.0),
      LinearTerm[LinearTerm("x2", 1.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Mixing consistent and inconsistent equations")
    add_variable!(rec, [x1, x2, x3, x4])
    add_equation!(rec, [eq1, eq2, eq3, eq4, eq5])

    res = remove_redundant_equations!(deepcopy(rec))
    @test iserror(res)
    @test issetequal(keys(res.error.infos[INFO_CONFLICTING_EQ_DEPS]), [4])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][4], [1, 3, 4])

    res = presolve(rec)
    @test iserror(res)
    @test issetequal(keys(res.error.infos[INFO_CONFLICTING_EQ_DEPS]), [4])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][4], [1, 3, 4])
  end

  @testset "Several conflicts" begin
    x1 = MeasuredVariable(name = "x1", value = 10.0, uncertainty = 1.0)
    x2 = UnmeasuredVariable(name = "x2")
    x3 = MeasuredVariable(name = "x3", value = 3.0, uncertainty = 0.3)
    x4 = MeasuredVariable(name = "x4", value = 5.0, uncertainty = 0.5)

    # x1 + x2 = -1
    eq1 = Equation(
      ConstantTerm(1.0),
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x2", 1.0)],
      BilinearTerm[],
    )

    # x2 + x3 = -1
    eq2 = Equation(
      ConstantTerm(1.0),
      LinearTerm[LinearTerm("x2", 1.0), LinearTerm("x3", 1.0)],
      BilinearTerm[],
    )

    # x1 - x3 = -1
    eq3 = Equation(
      ConstantTerm(1.0),
      LinearTerm[LinearTerm("x1", 1.0), LinearTerm("x3", -1.0)],
      BilinearTerm[],
    )

    # 3*x1 + 5*x2 + 2*x3 = -4
    eq4 = Equation(
      ConstantTerm(4.0),
      LinearTerm[LinearTerm("x1", 3.0), LinearTerm("x2", 5.0), LinearTerm("x3", 2.0)],
      BilinearTerm[],
    )

    rec = ReconciliationProblem(name = "Several conflicts")
    add_variable!(rec, [x1, x2, x3, x4])
    add_equation!(rec, [eq1, eq2, eq3, eq4])

    res = remove_redundant_equations!(deepcopy(rec))
    @test iserror(res)
    @test issetequal(keys(res.error.infos[INFO_CONFLICTING_EQ_DEPS]), [3, 4])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][3], [1, 2, 3])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][4], [1, 2, 4])

    res = presolve(rec)
    @test iserror(res)
    @test issetequal(keys(res.error.infos[INFO_CONFLICTING_EQ_DEPS]), [3, 4])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][3], [1, 2, 3])
    @test issetequal(res.error.infos[INFO_CONFLICTING_EQ_DEPS][4], [1, 2, 4])
  end
end

@testset "Unobservable variables" begin
  @testset "Check solving of unobservable variables" begin
    x1 = UnmeasuredVariable(; name = "x1", lb = -Inf)
    y1 = MeasuredVariable(name = "y1", value = 5.0, uncertainty = 0.5)
    z1 = FixedVariable(; name = "z1", value = 0.0)

    # Equation in the form of `2 + x1 + z1 * y1 = 0`.
    eq = parse_equation("2 + x1 + z1 * y1 = 0")
    problem = ReconciliationProblem()
    add_variable!(problem, [x1, y1, z1])
    add_equation!(problem, eq)

    result = unwrap(presolve(problem))
    @test result.replacements == dictionary(["z1" => 0.0, "x1" => -2.0])
    @test only(keys(result.unconstrained)) == "y1"

    sol = solve(problem)
    @test sol.status == FEASIBLE
    # Can only solve for 2 out of 3 variables.
    @test sol.score.solved_variables_score == 1
    @test sol.score.solved_equations_score == 1
  end

  @testset "Check identification of unobservable variables" begin
    @testset "Problem without equations" begin
      rec = ReconciliationProblem()
      add_variable!(rec, UnmeasuredVariable(; name = "x"))
      add_variable!(rec, UnmeasuredVariable(; name = "y"))
      prerec = unwrap(presolve(rec))

      @test isempty(prerec.pre.equations)
      @test isempty(prerec.pre.variables)
      @test isempty(prerec.variables_kept)
      @test issetequal(keys(prerec.unobservable), ["x", "y"])
    end

    @testset "Problem without equations after presolve" begin
      rec = ReconciliationProblem()
      add_variable!(rec, UnmeasuredVariable(; name = "x"))
      add_variable!(rec, UnmeasuredVariable(; name = "y"))
      add_variable!(rec, UnmeasuredVariable(; name = "z"))
      add_equation!(rec, Equation(ConstantTerm(-2.0), [LinearTerm("x", 1.0)], BilinearTerm[]))
      prerec = unwrap(presolve(rec))

      @test isempty(prerec.pre.equations)
      @test isempty(prerec.pre.variables)
      @test isempty(prerec.variables_kept)
      @test issetequal(keys(prerec.unobservable), ["y", "z"])
    end
  end
end

@testset "Augmented matrix" begin
  y1 = UnmeasuredVariable(; name = "y1")
  y2 = UnmeasuredVariable(; name = "y2")
  eq1 = parse_equation("2 * y1  + 3 * y2 = 5")
  eq2 = parse_equation("4 * y1 - y2 = 1")
  rec = ReconciliationProblem()
  add_variable!(rec, [y1, y2])
  add_equation!(rec, [eq1, eq2])

  # Check that the augmented matrix is recovered.
  A = [2 3; 4 -1]
  b = [5, 1]
  @test augmented_matrix(rec) == hcat(A, b)

  # Check correct solution is obtained.
  y1sol = 4 // 7
  y2sol = 9 // 7
  @test 2 * y1sol + 3 * y2sol == 5
  @test 4 * y1sol - y2sol == 1
  @test A \ b ≈ [y1sol, y2sol]
  sol = solve(rec)
  sol.variables_value["y1"] ≈ y1sol
  sol.variables_value["y2"] ≈ y2sol

  # The method is only defined for linear systems.
  eq3 = parse_equation("y1 * y2 = 1")
  add_equation!(rec, eq3)
  @test_throws ArgumentError augmented_matrix(rec)
end

@testset "independent_cols!" begin
  @testset "dense: two independent out of three columns" begin
    # Columns 1 and 2 are linearly dependent (col2 = 5 * col1).
    A = [1.0 5.0 0.0; 0.0 0.0 2.0; 1.0 5.0 0.0]
    cols, _ = independent_cols!(A, RREFConfig(atol = 1e-7))
    @test length(cols) == 2
    @test rank(A[:, cols]) == 2
  end

  @testset "dense: full rank" begin
    A = [1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0]
    cols, _ = independent_cols!(A, RREFConfig(atol = 1e-7))
    @test sort(cols) == [1, 2, 3]
  end

  @testset "sparse: two independent out of three columns" begin
    # Regression test: without applying QR.pcol the returned indices referred to
    # the permuted ordering rather than the original column indices of A.
    A = sparse([1.0 5.0 0.0; 0.0 0.0 2.0; 1.0 5.0 0.0])
    cols, _ = independent_cols!(A, SparseRREFConfig(atol = 1e-7))
    @test length(cols) == 2
    @test all(1 .<= cols .<= size(A, 2))
    @test rank(Matrix(A[:, cols])) == 2
  end

  @testset "sparse: full rank" begin
    A = sparse([1.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0])
    cols, _ = independent_cols!(A, SparseRREFConfig(atol = 1e-7))
    @test sort(cols) == [1, 2, 3]
  end

  @testset "sparse: all columns identical" begin
    A = sparse([1.0 1.0 1.0; 2.0 2.0 2.0])
    cols, _ = independent_cols!(A, SparseRREFConfig(atol = 1e-7))
    @test length(cols) == 1
  end
end

@testset "Test find_unobservable_variables_in_constraints with and without running presolve first with STAN models" begin
  for m in models
    @info "Running model: $(basename(m))"
    # Convert data to reconciliation problem.
    trace_filename = joinpath(m, "data", "trace.txt")
    if !isfile(trace_filename)
      @warn "Skipping model $(basename(m)). Expected trace file under `/data/trace.txt`."
      continue
    end
    strace = parse(STANTrace, m)
    rec = convert(ReconciliationProblem, strace)
    sol_stan = convert(ReconciliationSolution, strace)

    if isequal(sol_stan.status, INFEASIBLE) || isequal(sol_stan.status, CONTRADICTORY) ||
       isequal(get(sol_stan.infos, INFO_STAN_STATUS, nothing), "INFEASIBLE") ||
       isequal(get(sol_stan.infos, INFO_STAN_STATUS, nothing), "CONTRADICTORY")
      @info "Skipping model $(basename(m)) as it is infeasibile."
      continue
    end

    @time undet_vars = keys(find_unobservable_variables_in_constraints(sol_stan))
    all_vars = keys(rec.variables)
    sol_stan_vars = keys(sol_stan.variables_value)
    sol_stan_free = keys(sol_stan.unobservable_variables)

    # Check that no variable assigned by STAN is considered ambiguous without presolve
    if basename(m) in []
      test_broken = true
      @warn "Test that no variables assigned by STAN are considered free without presolve for model $(basename(m)) set to 'broken'."
    else
      test_broken = false
    end
    @test intersect(sol_stan_vars, undet_vars) == Indices{String}() broken = test_broken
    # Check that all variables not assigned by STAN are considered ambiguous without presolve
    if basename(m) in []
      test_broken = true
      @warn "Test that all variables not assigned by STAN are considered free without presolve for model $(basename(m)) set to 'broken'."
    else
      test_broken = false
    end
    @test issubset(sol_stan_free, undet_vars) broken = test_broken


    @time undet_vars_with_presolve = keys(find_unobservable_variables_in_constraints(sol_stan, run_presolve_first = true))

    # Check that no variable assigned by STAN is considered ambiguous using presolve
    if basename(m) in []
      test_broken = true
      @warn "Test that no variables assigned by STAN are considered free using presolve for model $(basename(m)) set to 'broken'."
    else
      test_broken = false
    end
    @test intersect(sol_stan_vars, undet_vars_with_presolve) == Indices{String}() broken = test_broken
    # Check that all variables not assigned by STAN are considered ambiguous using presolve
    if basename(m) in []
      test_broken = true
      @warn "Test that all variables not assigned by STAN are considered free using presolve for model $(basename(m)) set to 'broken'."
    else
      test_broken = false
    end
    @test issubset(sol_stan_free, undet_vars_with_presolve) broken = test_broken


    # Check that find_unobservable_variables_in_constraints with and without running presolve first identify the same variables as ambiguous
    if basename(m) in []
      test_broken = true
      @warn "Test that unobservable variable identification is independent of presolving before for model $(basename(m)) set to 'broken'."
    else
      test_broken = false
    end
    @test issetequal(undet_vars, undet_vars_with_presolve) broken = test_broken
  end
end

@testset "Identification of unobservable variables" begin
  @testset "Two unknown, one equation" begin
    x = FixedVariable(; name = "x", value = 2.0)
    y = UnmeasuredVariable(; name = "y")
    z = UnmeasuredVariable(; name = "z")
    eq = parse_equation("x + y = z")
    problem = ReconciliationProblem()
    add_variable!(problem, [x, y, z])
    add_equation!(problem, eq)

    free_vars = find_unobservable_variables_in_constraints(problem, Dictionary(["x", "y", "z"], [2.0, 1.0, 3.0]))
    @test issetequal(keys(free_vars), ["y", "z"])
  end

  @testset "Unique definition due to variable bounds" begin
    i = FixedVariable(; name = "i", value = 1.0)
    o1 = UnmeasuredVariable(; name = "o1")
    o2 = UnmeasuredVariable(; name = "o2")
    o3 = UnmeasuredVariable(; name = "o3")
    tc1 = FixedVariable(; name = "tc1", value = 1.0, tc = true)
    tc2 = UnmeasuredVariable(; name = "tc2", tc = true)
    tc3 = UnmeasuredVariable(; name = "tc3", tc = true)
    eq1 = parse_equation("o1 = tc1 * i")
    eq2 = parse_equation("o2 = tc2 * i")
    eq3 = parse_equation("o3 = tc3 * i")
    eq4 = parse_equation("tc1 + tc2 + tc3 = 1.0")
    problem = ReconciliationProblem()
    add_variable!(problem, [i, o1, o2, o3, tc1, tc2, tc3])
    add_equation!(problem, [eq1, eq2, eq3, eq4])
    sol = solve(problem)

    expected_values = Dict("i" => 1.0, "o1" => 1.0, "o2" => 0.0, "o3" => 0.0, "tc1" => 1.0, "tc2" => 0.0, "tc3" => 0.0)
    for var in keys(expected_values)
      @test relaxed_isapprox(sol.variables_value[var], expected_values[var])
    end

    @test isempty(sol.unobservable_variables)
  end

  @testset "Unobservable variables identification for numerically sensitive problem" begin
    filepath = joinpath(pkgdir(OpenSEFA), "test", "examples", "Cencic_2026_1-model-v1.json")
    @assert isfile(filepath)
    problem = load_model(filepath)

    filepath2 = joinpath(pkgdir(OpenSEFA), "test", "examples", "Cencic_2026_1-model-v2.json")
    @assert isfile(filepath2)
    problem2 = load_model(filepath2)

    sol_default = solve(problem)
    sol_jump = solve(problem, JuMPSolver())

    sol2_default = solve(problem2)
    sol2_jump = solve(problem2, JuMPSolver())

    @test issetequal(keys(sol_default.unobservable_variables), keys(sol_jump.unobservable_variables))
    @test issetequal(keys(sol2_default.unobservable_variables), keys(sol2_jump.unobservable_variables))
    @test issetequal(keys(sol_default.unobservable_variables), keys(sol2_default.unobservable_variables))
  end

  @testset "Seven flows variant" begin
    # Variant of the seven-flows example where m6 and m7 are unobservable.
    # The system has 4 equations in 8 variables (m1..m7, tc34):
    #
    #   eq1:  m1 + m2 + m4 - m3 = 0   (mass balance at node 1)
    #   eq2:  m3 - m4 - m5 = 0        (mass balance at node 2)
    #   eq3:  m5 - m6 - m7 = 0        (mass balance at node 3)
    #   eq4:  m4 - m3 * tc34 = 0      (transfer coefficient relation)
    #
    # m1, m2, m5 are fixed and m4 is determined by eq4 after presolving tc34.
    # eq3 constrains only the sum m6 + m7 = m5, leaving one degree of freedom,
    # so m6 and m7 are individually unobservable.
    m1 = FixedVariable(; name = "m1", value = 100.0)
    m2 = FixedVariable(; name = "m2", value = 50.0)
    m3 = MeasuredVariable(; name = "m3", value = 300.0, uncertainty = 30.0)
    m4 = UnmeasuredVariable(; name = "m4")
    m5 = FixedVariable(; name = "m5", value = 150.0)
    m6 = UnmeasuredVariable(; name = "m6")
    m7 = UnmeasuredVariable(; name = "m7")
    tc34 = MeasuredVariable(; name = "tc34", value = 0.5, uncertainty = 0.5)

    eq1 = parse_equation("m1 + m2 + m4 - m3 = 0")
    eq2 = parse_equation("m3 - m4 - m5 = 0")
    eq3 = parse_equation("m5 - m6 - m7 = 0")
    eq4 = parse_equation("m4 - m3 * tc34 = 0")

    problem = ReconciliationProblem()
    add_variable!(problem, [m1, m2, m3, m4, m5, m6, m7, tc34])
    add_equation!(problem, [eq1, eq2, eq3, eq4])

    sol = solve(problem, JuMPSolver(; tol = 1e-10, max_iter = 10000))
    free_vars = collect(keys(sol.unobservable_variables))
    num_free = sol.infos[INFO_N_UNOBSERVABLE]
    @test num_free == 2
    @test Set(free_vars) == Set([m6.name, m7.name])
  end
end
