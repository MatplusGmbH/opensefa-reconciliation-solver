using Test
using OpenSEFA.STANModule
using OpenSEFA.ConstraintsModule
using OpenSEFA.ReconciliationModule
import OpenSEFA as OSEFA
using ResultTypes: iserror, unwrap

# Filepath for each STAN model available under the `models` folder.
# Returns empty when MODELS_PATH doesn't exist; model-dependent tests skip gracefully.
models = all_models()

@testset "Loading STAN CSV data." begin
  for m in models
    # Check that all available model data can be loaded into Julia data structures.
    # The STAN .csv files should be located under `/model/data` for each model.
    data_path = joinpath(m, "data")
    if !isdir(data_path)
      @warn ArgumentError("CSV data not available under $data_path.")
      continue
    end
    data = safe_load_stan_data(m)
    if iserror(data)
      @info "CSV data for model $m not found."
      continue
    end
    data = unwrap(data)
    @test data isa STANData
  end
end

@testset "Parsing STAN Trace Calculation data." begin
  for m in models
    trace_filename = joinpath(m, "data", "trace.txt")
    if !isfile(trace_filename)
      @warn "Skipping model $(basename(m)). Expected trace file under `/data/trace.txt`."
      continue
    end
    data = parse(STANTrace, m)
    @test data isa STANTrace
  end
end

@testset "Parsing result from Kelly solver." begin
  kelly_models = [
    "4_PVC Austria (1950-1994) plus 1 Layer",
    "4_PVC Austria (1950-1994) plus 2 Layers",
    "Status Quo",
  ]
  models_dir = joinpath(pkgdir(OSEFA), "models")
  if !isdir(models_dir)
    @warn "Skipping Kelly solver test: models directory not found at $models_dir."
  else
    for m in kelly_models
      @info "Running model : $(basename(m))"
      filepath = joinpath(models_dir, m)
      trace = parse(STANTrace, filepath)
      rec = convert(ReconciliationProblem, trace)

      @test !ismissing(trace.kelly_results)
      @test length(trace.kelly_results.variables_value) + length(trace.kelly_results.unassigned_variables) ==
            length(rec.variables)
    end
  end
end

@testset "Parsing STAN in German." begin
  @testset "Loading STAN CSV data." begin
    m = joinpath(pkgdir(OSEFA), "test", "examples", "STAN_German")
    data_path = joinpath(m, "data")
    if !isdir(data_path)
      @warn ArgumentError("CSV data not available under $data_path.")
    end
    data = safe_load_stan_data(m)
    if iserror(data)
      @info "CSV data for model $m not found."
    else
      data = unwrap(data)
      @test data isa STANData
    end
  end

  @testset "Parsing STAN Trace Calculation data." begin
    m = joinpath(pkgdir(OSEFA), "test", "examples", "STAN_German")
    trace_filename = joinpath(m, "data", "trace.txt")
    if !isfile(trace_filename)
      @warn "Skipping model $(basename(m)). Expected trace file under `/data/trace.txt`."
    end
    data = parse(STANTrace, m)
    @test data isa STANTrace
    @test length(data.initial_variables) == 13
    @test data.initial_variables[9].name == "TC:18"
    @test data.initial_variables[9].value == 0.1
    @test data.initial_variables[9].uncertainty == 0.1
    @test data.initial_variables[9].tc
    @test length(data.equations) == 9
    @test data.equations[1] == Equation(ConstantTerm(), [LinearTerm("BS:3", 0.5)], BilinearTerm[])
  end
end
