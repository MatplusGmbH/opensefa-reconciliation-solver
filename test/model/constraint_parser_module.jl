
using Test
using OpenSEFA.ConstraintsModule
using OpenSEFA.ConstraintParserModule

@testset "Simple assignments" begin
  # Test R = 100.0 → R - 100.0 = 0
  result = parse_constraint("R = 100.0")
  @test isa(result, Equation)
  @test result.constant_term.value == -100.0  # R - 100 = 0, so constant is -100
  @test only(result.linear_terms) == LinearTerm("R", 1.0)
  @test isempty(result.bilinear_terms)

  # Test E = 400 → E - 400 = 0
  result = parse_constraint("E = 400")
  @test isa(result, Equation)
  @test result.constant_term == ConstantTerm(-400.0)
  @test only(result.linear_terms) == LinearTerm("E", 1.0)
  @test isempty(result.bilinear_terms)
end

@testset "Bilinear assignments" begin
  # Test F1_mass_N1 = R*x1 → F1_mass_N1 - R*x1 = 0
  result = parse_constraint("F1_mass_N1 = R*x1")
  @test isa(result, Equation)
  @test result.constant_term == ConstantTerm(0.0)
  @test only(result.linear_terms) == LinearTerm("F1_mass_N1", 1.0)
  # Negated because moved from RHS.
  @test only(result.bilinear_terms) == BilinearTerm("R", "x1", -1.0)

  # Test division by numbers
  result = parse_constraint("F1_mass_N1 = R*x1 + 10/3")
  @test result.constant_term == ConstantTerm(-10 / 3)
  result = parse_constraint("4 = 6 * x1 / 3")
  @test result.linear_terms[1] == LinearTerm("x1", -2.0)
  result = parse_constraint("4 = 6 * x1 / 3 / 2")
  @test result.linear_terms[1] == LinearTerm("x1", -1.0)

  # Test a = -b + c
  result = parse_constraint("a = -b + c")
  @test isa(result, Equation)
  @test result.constant_term == ConstantTerm(0.0)
  @test result.linear_terms == [
    LinearTerm("a", 1.0),
    LinearTerm("b", 1.0),
    LinearTerm("c", -1.0),
  ]
  @test isempty(result.bilinear_terms)
end

@testset "Linear constraints" begin
  # Test backward compatibility: x + y = 0
  result = parse_constraint("x + y = 0")
  @test isa(result, Equation)
  @test result.constant_term == ConstantTerm(0.0)
  @test result.linear_terms == [LinearTerm("x", 1.0), LinearTerm("y", 1.0)]
  @test isempty(result.bilinear_terms)

  # Test x1 + x2 + x3 = 1 → x1 + x2 + x3 - 1 = 0
  result = parse_constraint("x1 + x2 + x3 = 1")
  @test isa(result, Equation)
  # -1 from moving RHS to LHS
  @test result.constant_term.value == -1.0
  @test result.linear_terms == [
    LinearTerm("x1", 1.0)
    LinearTerm("x2", 1.0)
    LinearTerm("x3", 1.0)
  ]
  @test isempty(result.bilinear_terms)
end

@testset "Edge cases and error handling" begin
  # Invalid constraint format.
  @test_throws ArgumentError parse_constraint("invalid constraint")

  # Unrecognized term structure.
  # Too many components.
  @test_throws ArgumentError parse_constraint("x = a*b*c*d")

  # Unsupported term structure.
  # Dividing by variable.
  @test_throws ArgumentError parse_constraint("x = a/b")

  # Invalid bound constraint.
  # Non-numeric bounds.
  @test_throws ArgumentError parse_constraint("x <= y <= z")
end

@testset "Bound constraints" begin
  test_cases = [
    # Double bounds - ASCII operators
    "0 <= x <= 1" => BoundConstraint("x", 0.0, 1.0),
    "0<x<1" => BoundConstraint("x", 0.0, 1.0),
    "0 <= x < 1" => BoundConstraint("x", 0.0, 1.0),
    "0 < x <= 1" => BoundConstraint("x", 0.0, 1.0),
    "-5 <= y <= 5" => BoundConstraint("y", -5.0, 5.0),
    "0.1 <= z <= 0.9" => BoundConstraint("z", 0.1, 0.9),
    "0 < Cu_steel <= 0.15" => BoundConstraint("Cu_steel", 0.0, 0.15),

    # Double bounds - Unicode operators
    "0 ≤ x ≤ 1" => BoundConstraint("x", 0.0, 1.0),
    "0 ≤ temp ≤ 100" => BoundConstraint("temp", 0.0, 100.0),

    # Single bounds - variable op number format
    "x >= 0" => BoundConstraint("x", 0.0, missing),
    "x > 0" => BoundConstraint("x", 0.0, missing),
    "x <= 100" => BoundConstraint("x", missing, 100.0),
    "x < 100" => BoundConstraint("x", missing, 100.0),
    "y >= -10" => BoundConstraint("y", -10.0, missing),
    "z <= 3.14" => BoundConstraint("z", missing, 3.14),

    # Single bounds - Unicode operators
    "x ≥ 0" => BoundConstraint("x", 0.0, missing),
    "x ≤ 100" => BoundConstraint("x", missing, 100.0),

    # Single bounds - number op variable format
    "0 <= x" => BoundConstraint("x", 0.0, missing),
    "0 < x" => BoundConstraint("x", 0.0, missing),
    "100 >= x" => BoundConstraint("x", missing, 100.0),
    "100 > x" => BoundConstraint("x", missing, 100.0),
    "-10 <= y" => BoundConstraint("y", -10.0, missing),
    "3.14 >= z" => BoundConstraint("z", missing, 3.14),

    # Single bounds - Unicode number op variable format
    "0 ≤ x" => BoundConstraint("x", 0.0, missing),
    "100 ≥ x" => BoundConstraint("x", missing, 100.0),

    # Complex variable names
    "steel_density >= 7.85" => BoundConstraint("steel_density", 7.85, missing),
    "F1_mass_N1 > 0" => BoundConstraint("F1_mass_N1", 0.0, missing),
    "FV:726 <= 1000" => BoundConstraint("FV:726", missing, 1000.0),
    "hourly_consumption <= 1000.5" => BoundConstraint("hourly_consumption", missing, 1000.5),

    # Edge cases with decimals
    "pressure >= 101.325" => BoundConstraint("pressure", 101.325, missing),
    "efficiency <= 0.95" => BoundConstraint("efficiency", missing, 0.95),
    "temperature >= -273.15" => BoundConstraint("temperature", -273.15, missing),

    # Mixed Unicode and ASCII
    "0 ≤ x <= 1" => BoundConstraint("x", 0.0, 1.0),
    "0 <= x ≤ 1" => BoundConstraint("x", 0.0, 1.0),
  ]
  for (expr, expected) in test_cases
    @test parse_constraint(expr) == expected
  end
end

@testset "More complex constraints" begin
  # Test scrap_cost = F1_mass_N1 * C_Mass_F1 + F2_mass_N1 * C_Mass_F2 + F3_mass_N1 * C_Mass_F3
  # → scrap_cost - F1_mass_N1 * C_Mass_F1 - F2_mass_N1 * C_Mass_F2 - F3_mass_N1 * C_Mass_F3 = 0
  result = parse_constraint("scrap_cost = F1_mass_N1 * C_Mass_F1 + F2_mass_N1 * C_Mass_F2 + F3_mass_N1 * C_Mass_F3")
  @test isa(result, Equation)
  @test result.constant_term == ConstantTerm(0.0)
  @test only(result.linear_terms) == LinearTerm("scrap_cost", 1.0)
  # All bilinear terms should have coefficient -1.0 (moved from RHS).
  @test result.bilinear_terms == [
    BilinearTerm("F1_mass_N1", "C_Mass_F1", -1.0)
    BilinearTerm("F2_mass_N1", "C_Mass_F2", -1.0)
    BilinearTerm("F3_mass_N1", "C_Mass_F3", -1.0)
  ]

  # Test total_cost = scrap_cost + energy_cost
  # → total_cost - scrap_cost - energy_cost = 0
  result = parse_constraint("total_cost = scrap_cost + energy_cost")
  @test isa(result, Equation)
  @test result.constant_term == ConstantTerm(0.0)
  @test result.linear_terms == [
    LinearTerm("total_cost", 1.0)
    LinearTerm("scrap_cost", -1.0)
    LinearTerm("energy_cost", -1.0)
  ]
  @test isempty(result.bilinear_terms)
end

@testset "Constraint with numerical values" begin
  # Check that the numeric value can be defined both on the left and on the right of the term.
  eq = "F7_energy_N1 - 1.5 * F9_energy_N1 = 0"
  cons = parse_constraint(eq)
  @test isa(cons, Equation)
  @test cons.constant_term == ConstantTerm(0.0)
  @test cons.linear_terms == [LinearTerm("F7_energy_N1", 1.0), LinearTerm("F9_energy_N1", -1.5)]
  @test isempty(cons.bilinear_terms)

  eq = "F7_energy_N1 - F9_energy_N1 * 1.5 = 0"
  cons = parse_constraint(eq)
  @test isa(cons, Equation)
  @test cons.constant_term == ConstantTerm(0.0)
  @test cons.linear_terms == [LinearTerm("F7_energy_N1", 1.0), LinearTerm("F9_energy_N1", -1.5)]
  @test isempty(cons.bilinear_terms)
end

@testset "Equations with zero factor terms" begin
  eq1 = parse_equation("3 * x + 2*x *y = 12")
  @test eq1.constant_term == ConstantTerm(-12.0)
  @test eq1.linear_terms == [LinearTerm("x", 3.0)]
  @test eq1.bilinear_terms == [BilinearTerm("x", "y", 2.0)]

  eq2 = parse_equation("0 * x + 2*x *y = 12")
  @test eq2.constant_term == ConstantTerm(-12.0)
  @test eq2.linear_terms == LinearTerm[]
  @test eq2.bilinear_terms == [BilinearTerm("x", "y", 2.0)]

  eq3 = parse_equation("3 * x + 0*x *y = 12")
  @test eq3.constant_term == ConstantTerm(-12.0)
  @test eq3.linear_terms == [LinearTerm("x", 3.0)]
  @test eq3.bilinear_terms == BilinearTerm[]
end
