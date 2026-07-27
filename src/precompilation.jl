using Suppressor: @suppress_out

# Precompilation workload for OpenSEFA.
# Uses only inline model definitions to avoid depending on external data files.
@setup_workload begin
  @compile_workload begin
    include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))

    rec = seven_flows_cencic_2012()

    @suppress_out begin
      solve(rec)
    end
  end
end
