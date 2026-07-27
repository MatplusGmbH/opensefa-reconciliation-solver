using Test
using OpenSEFA.ComparisonModule

@testset "Tolerances" begin
  rtol = get_rtol(Float64)
  ztol = get_ztol(Float64)
  atol = get_atol(Float64)
  @test rtol isa Float64 && rtol > 0
  @test ztol isa Float64 && ztol > 0
  @test atol isa Float64 && atol >= 0
end

@testset "isapproxzero" begin
  @test isapproxzero(0.0)
  @test isapproxzero(1e-16)
  @test !isapproxzero(1.0)
  @test isapproxzero(0.0; ztol = 1e-10)
  @test !isapproxzero(1e-9; ztol = 1e-10)
end

@testset "relaxed_isapprox" begin
  @test relaxed_isapprox(1.0, 1.0)
  @test relaxed_isapprox(1.0, 1.0 + 1e-8)
  @test !relaxed_isapprox(1.0, 2.0)
  @test relaxed_isapprox(0.0, 0.0)
  # Mixed types
  @test relaxed_isapprox(1.0, 1)
end

@testset "relaxed_leq and relaxed_geq" begin
  @test relaxed_leq(1.0, 2.0)
  @test relaxed_leq(1.0, 1.0)
  @test !relaxed_leq(2.0, 1.0)
  @test relaxed_geq(2.0, 1.0)
  @test relaxed_geq(1.0, 1.0)
  @test !relaxed_geq(1.0, 2.0)
end

@testset "tight_leq and tight_geq" begin
  @test tight_leq(1.0, 2.0)
  @test !tight_leq(2.0, 1.0)
  @test tight_geq(2.0, 1.0)
  @test !tight_geq(1.0, 2.0)
end

@testset "isapprox_dict" begin
  using Dictionaries
  d1 = Dictionary(["a", "b"], [1.0, 2.0])
  d2 = Dictionary(["a", "b"], [1.0, 2.0])
  @test isapprox_dict(d1, d2)
  d3 = Dictionary(["a", "b"], [1.0, 3.0])
  @test !isapprox_dict(d1, d3)
end
