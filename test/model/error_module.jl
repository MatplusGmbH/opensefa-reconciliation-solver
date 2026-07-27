using Test
using Dictionaries
using OpenSEFA.ErrorModule
using OpenSEFA.ErrorInfosModule

@testset "Error type hierarchy" begin
  @test ContradictionError <: SolverError
  @test InfeasibleEquationError <: ContradictionError
  @test IllposedProblemError <: ContradictionError
  @test InternalError <: SolverError
  @test FileIOError <: Exception
end

@testset "InfeasibleEquationError" begin
  err = InfeasibleEquationError("value out of range"; index=[1], residual=[0.0])
  @test err isa ContradictionError
  @test err isa SolverError
  @test_throws InfeasibleEquationError throw(err)
end

@testset "IllposedProblemError" begin
  err = IllposedProblemError("duplicate variable name 'x'")
  @test err isa ContradictionError
  @test_throws IllposedProblemError throw(err)
end

@testset "InternalError" begin
  err = InternalError("unexpected state")
  @test err isa SolverError
  @test_throws InternalError throw(err)
end

@testset "FileIOError" begin
  err = FileIOError("file not found: model.json")
  @test err isa Exception
  @test_throws FileIOError throw(err)
end

@testset "add_info!" begin
  info = Dictionary{String, Any}()
  add_info!(info, "key1", 42)
  @test info["key1"] == 42
  # Adding new key works
  add_info!(info, "key2", "hello")
  @test info["key2"] == "hello"
  # Adding to an existing key accumulates into a vector
  add_info!(info, "key1", 99)
  @test info["key1"] == [42, 99]
end
