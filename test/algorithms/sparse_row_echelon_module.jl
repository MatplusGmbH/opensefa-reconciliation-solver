using Test
using OpenSEFA.SparseRowEchelonModule
using SparseArrays

@testset "RREF with Pivots Tests" begin

  @testset "Simple dependence." begin
    A = [1 2 3; 2 4 5.0]
    expected_rref = [1.0 2.0 0.0; 0.0 0.0 1.0]
    expected_pivots = [1, 3]

    rref, pivots = sparse_rref_with_pivots(sparse(A))
    @test isapprox(rref, expected_rref, atol = 1e-8)
    @test pivots == expected_pivots
  end

  @testset "Mixed linear independence/dependence." begin
    B = [1 0 1 2 1; 0 1 1 1 1; 1 1 2 3 1.0]
    expected_rref = [
      1 0 1 2 0;
      0 1 1 1 0;
      0 0 0 0 1
    ]
    expected_pivots = [1, 2, 5]

    rref, pivots = sparse_rref_with_pivots(sparse(B))
    @test isapprox(rref, expected_rref, atol = 1e-8)
    @test pivots == expected_pivots
  end

  @testset "4D test case 1." begin
    A1 = [1 2 -1 -4; 2 3 -1 -11; -2 0 -3 22.0]
    expected_rref_1 = [
      1.0 0.0 0.0 -8.0;
      0.0 1.0 0.0 1.0;
      0.0 0.0 1.0 -2.0
    ]
    expected_pivots_1 = [1, 2, 3]

    rref_1, pivots_1 = sparse_rref_with_pivots(sparse(A1))
    @test isapprox(rref_1, expected_rref_1, atol = 1e-8)
    @test pivots_1 == expected_pivots_1
  end

  @testset "4D test case 2." begin
    A2 = [16 2 3 13; 5 11 10 8; 9 7 6 12; 4 14 15 1.0]
    expected_rref_2 = [
      1.0 0.0 0.0 1.0;
      0.0 1.0 0.0 3.0;
      0.0 0.0 1.0 -3.0;
      0.0 0.0 0.0 0.0
    ]
    expected_pivots_2 = [1, 2, 3]

    rref_2, pivots_2 = sparse_rref_with_pivots(sparse(A2))
    @test isapprox(rref_2, expected_rref_2, atol = 1e-8)
    @test pivots_2 == expected_pivots_2
  end

  @testset "4d test case 3." begin
    A3 = [1 2 0 3; 2 4 0 7.0]
    expected_rref_3 = [
      1.0 2.0 0.0 0.0;
      0.0 0.0 0.0 1.0
    ]
    expected_pivots_3 = [1, 4]

    rref_3, pivots_3 = sparse_rref_with_pivots(sparse(A3))
    @test isapprox(rref_3, expected_rref_3, atol = 1e-8)
    @test pivots_3 == expected_pivots_3
  end
end
