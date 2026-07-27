using Test
import OpenSEFA
using OpenSEFA.ConstraintsModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.PresolveModule: remove_fixed_variables!

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))

@testset "Manual model creation." begin
  rec = seven_flows_cencic_2012()
  @test rec isa ReconciliationProblem

  # Check removal of fixed variable.
  crec = deepcopy(rec)
  remove_fixed_variables!(crec)
  @test isempty(crec.fixed)
end
