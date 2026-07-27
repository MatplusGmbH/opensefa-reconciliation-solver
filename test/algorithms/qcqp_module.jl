using Test
import OpenSEFA
using OpenSEFA.QCQPModule
using OpenSEFA.ConstraintsModule
using OpenSEFA.ReconciliationModule

# Build the seven-flows problem for conversion tests.
include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))

@testset "QCQP conversion from ReconciliationProblem" begin
  rec = seven_flows_cencic_2012()
  qcqp = convert(QuadraticallyConstrainedQuadraticProgram, rec)
  @test qcqp isa QuadraticallyConstrainedQuadraticProgram
  @test QCQP === QuadraticallyConstrainedQuadraticProgram
  # One QCQP constraint per reconciliation equation.
  @test length(qcqp.constraints) == length(rec.equations)
end

@testset "QuadraticConstraint evaluation" begin
  rec = seven_flows_cencic_2012()
  qcqp = convert(QuadraticallyConstrainedQuadraticProgram, rec)
  n = length(qcqp.variable_index)
  u0 = zeros(n)
  # Evaluate at the initial point — result is a scalar.
  for c in qcqp.constraints
    val = evaluate_unsafe(c, u0)
    @test val isa Float64
  end
end

@testset "QuadraticConstraint gradient" begin
  rec = seven_flows_cencic_2012()
  qcqp = convert(QuadraticallyConstrainedQuadraticProgram, rec)
  n = length(qcqp.variable_index)
  u0 = zeros(n)
  grad = zeros(n)
  for c in qcqp.constraints
    fill!(grad, 0.0)
    evaluate_gradient_unsafe!(grad, c, u0)
    @test length(grad) == n
  end
end
