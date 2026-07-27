"""
Methods to parse STAN model data and calculations output.

STAN is available at https://www.stan2web.net/
"""
module STANModule

import CSV
import JSON3
using JSON3: write
using DataFrames: DataFrame
using Dictionaries
using ResultTypes: Result, @try, iserror, unwrap
using StringEncodings

import ..OpenSEFA: OpenSEFA
using ..ConstraintsModule
using ..ErrorModule
using ..ErrorInfosModule
using ..ConstraintParserModule
using ..ComparisonModule

export
  load_stan_data,
  safe_load_stan_data,
  all_models,
  write_json,
  MODELS_PATH,
  parse_kelly_trace,
  parse_variable!,
  parse_equation!,
  STANData,
  STANResult,
  STANTrace,
  STANInfos,
  add_STAN_infos!

"Path for each available model."
const MODELS_PATH = joinpath(pkgdir(OpenSEFA), "models")

"These models are work-in-progress, hence aren't tested yet."
IGNORED_MODELS = String[
  "4_PVC Austria (1950-1994) plus 1 Layer",
  "4_PVC Austria (1950-1994) plus 2 Layers",
  "Status Quo",
]

"""
Return the directory path of all available STAN models.
"""
function all_models(; path::String = MODELS_PATH, ignored_models::Vector{String} = IGNORED_MODELS)::Vector{String}
  model_path = String[]
  for (root, dirs, files) in walkdir(path)
    for dir in dirs
      # Search if the model file exists.
      candidate = dir * ".zmfa"
      for (_, _, f) in walkdir(joinpath(root, dir))
        if candidate ∈ f
          push!(model_path, joinpath(root, dir))
          break
        end
      end
    end
  end
  setdiff(model_path, joinpath.(MODELS_PATH, ignored_models))
end

# ============================
# STAN Data explorer parsing
# ============================

"Calculation algorithm in STAN used to solve the problem."
@enum STANCalculationMethod CENCIC_2012 IAL_IAMPL_2013

"""
Struct which holds data in CSV format produced by STAN software.
"""
struct STANData
  "Model's name."
  name::String

  "Path to the model data."
  path::String

  "Method used to solve the data reconciliation problem for the calculated values."
  method::STANCalculationMethod

  """
  DataFrame containing information about flows between processes.
  Columns:
  - Short symbol: Abbreviated name of the flow.
  - Name: Full name of the flow.
  - Remarks: Additional notes or comments.
  - From: Source of the flow.
  - To: Destination of the flow.
  """
  flows::DataFrame

  """
  DataFrame containing the measured and calculated values for flows.
  Columns:
  - Period, Layer: Time and structural layers.
  - Flow, Flow name: Identifiers for the flow.
  - Mass flow, ± Mass flow: Measured mass flow and its uncertainty.
  - Mass fraction, ± Mass fraction: Measured mass fraction and its uncertainty.
  - Mass flow (calculated), ± Mass flow (calculated): Reconciled mass flow and its uncertainty.
  - Volume flow, ± Volume flow: Measured volume flow and its uncertainty.
  - Mass concentration, ± Mass concentration: Measured mass concentration and its uncertainty.
  - Volume flow (calculated), ± Volume flow (calculated): Reconciled volume flow and its uncertainty.
  - Mass concentration (calculated), ± Mass concentration (calculated): Reconciled mass concentration and its uncertainty.
  - Literature reference, Remarks: References and additional notes.
  """
  flowvalues::DataFrame

  """
  DataFrame describing the processes involved in the system.
  Columns:
  - Short symbol: Abbreviated name of the process.
  - Name: Full name of the process.
  - Remarks: Additional notes or comments.
  - Stock: Stock information related to the process.
  - TCs: Transfer coefficients for the process.
  - Sub system: Whether the process describes a subsystem (boolean).
  """
  processes::DataFrame

  """
  DataFrame describing changes in stock values over time.
  Columns:
  - Period, Layer: Time and structural layers.
  - Process, Process name: Identifiers for the process.
  - Mass flow, ± Mass flow: Measured mass flow and its uncertainty for stock changes.
  - Mass fraction, ± Mass fraction: Measured mass fraction and its uncertainty for stock changes.
  - Mass flow (calculated), ± Mass flow (calculated): Reconciled mass flow and its uncertainty for stock changes.
  - Volume flow, ± Volume flow: Measured volume flow and its uncertainty for stock changes.
  - Mass concentration, ± Mass concentration: Measured mass concentration and its uncertainty for stock changes.
  - Volume flow (calculated), ± Volume flow (calculated): Reconciled volume flow and its uncertainty for stock changes.
  - Mass concentration (calculated), ± Mass concentration (calculated): Reconciled mass concentration and its uncertainty for stock changes.
  """
  stockdelta::DataFrame

  """
  DataFrame containing stock values for processes.
  Columns:
  - Period, Layer: Time and structural layers.
  - Process, Process name: Identifiers for the process.
  - Mass, ± Mass: Measured mass and its uncertainty.
  - Mass fraction, ± Mass fraction: Measured mass fraction and its uncertainty.
  - Stock (calculated), ± Stock (calculated): Reconciled stock value and its uncertainty.
  - Volume, ± Volume: Measured volume and its uncertainty.
  - Mass concentration, ± Mass concentration: Measured mass concentration and its uncertainty.
  - Stock volume (calculated), ± Stock volume (calculated): Reconciled stock volume and its uncertainty.
  - Mass concentration (calculated), ± Mass concentration (calculated): Reconciled mass concentration and its uncertainty.
  - Literature reference, Remarks: References and additional notes.
  """
  stockvalues::DataFrame

  """
  DataFrame describing system-wide settings and scaling information.
  Columns:
  - Name: Name of the system.
  - Scaling unit: Unit used for scaling.
  - Scaling factor: Factor used for scaling.
  - Description: Description of the system.
  """
  system::DataFrame

  """
  DataFrame describing transfer coefficients for processes.
  Columns:
  - Period, Layer: Time and structural layers.
  - Process, Process name: Identifiers for the process.
  - In->Out: Transfer pathway from input to output.
  - TC, ± TC: Transfer coefficient and its uncertainty.
  - TC (calculated), ± TC (calculated): Reconciled transfer coefficient and its uncertainty.
  - Literature reference, Remarks: References and additional notes.
  """
  transfer_coefficients::DataFrame
end

"""
Load the STAN data from `.csv` files located in `path`.
Such files are assumed to be under the `path/data` subfolder.
"""
function safe_load_stan_data(path::String)::Result{STANData, FileIOError}
  data_path = joinpath(path, "data")
  if !isdir(data_path)
    throw(ArgumentError("CSV data not available under $data_path."))
  end

  function read_df(file)::Result{DataFrame, FileIOError}
    filepath = joinpath(data_path, file * ".csv")
    if !isfile(filepath)
      return FileIOError("File $filepath not found.")
    end
    CSV.read(open(filepath, enc"UTF-16"), DataFrame)
  end

  """
  Try to read a DataFrame with English filename first, then German filename if English fails.
  """
  function read_df_bilingual(en_name::String, de_name::String)::Result{DataFrame, FileIOError}
    result = read_df(en_name)
    iserror(result) ? read_df(de_name) : result
  end

  name = basename(path)
  flows = @try read_df_bilingual("Flows", "Flüsse")
  flowvalues = @try read_df_bilingual("Flowvalues", "Flusswerte")
  processes = @try read_df_bilingual("Processes", "Prozesse")
  stockdelta = @try read_df_bilingual("Stockdelta", "Lageränderung")
  stockvalues = @try read_df_bilingual("Stockvalues", "Lagerwerte")
  system = @try read_df("System")
  transfer_coefficients = @try read_df_bilingual("Transfercoefficients", "Transferkoeffizienten")

  # Assume unless otherwise specified.
  method = CENCIC_2012

  STANData(
    name,
    data_path,
    method,
    flows,
    flowvalues,
    processes,
    stockdelta,
    stockvalues,
    system,
    transfer_coefficients,
  )
end
function load_stan_data(path::String)::STANData
  @try safe_load_stan_data(path)
end

"Serialize the data into JSON format."
function write_json(data::STANData)
  json_path = joinpath(dirname(data.path), data.name * ".json")
  open(json_path, "w") do io
    JSON3.write(io, data)
  end
  @info "Data stored in $json_path."
end


# ====================
# STAN Infos parsing
# ====================

"""
Struct which holds solving data from STAN
"""
struct STANInfos
  "Number of equations"
  neqs::Union{Int64, Nothing}
  "Number of variables"
  nvars::Union{Int64, Nothing}
  "Status"
  stan_status::Union{String, Nothing}
  "Modelling time (s)"
  modelling_time::Union{Float64, Nothing}
  "Computation time (s)"
  computation_time::Union{Float64, Nothing}
  "STAN version"
  version::Union{String, Nothing}
  "Machine information"
  machine_infos::Dictionary{String, String}
end

@enum InfosSection NEqsSection NVarsSecion StatusSection ModellingTimeSection ComputationTimeSection VersionSection MachineInfosSection

"Parse solving data from STAN"
function Base.parse(::Type{STANInfos}, path::String)::STANInfos
  filename = joinpath(path, "data", "solving_data.txt")
  name = basename(path)
  if !isfile(filename)
    throw(ArgumentError("Expected solving data file under `/data/solving_data.txt`."))
  end

  neqs = nothing
  nvars = nothing
  stan_status = nothing
  modelling_time = nothing
  computation_time = nothing
  version = nothing
  machine_infos = Dictionary{String, String}()

  fileresult = open(filename, "r") do file
    sec = nothing

    for line in eachline(file)
      if isempty(line)
        continue
      elseif occursin("Number of equations", line)
        sec = NEqsSection
      elseif occursin("Number of variables", line)
        sec = NVarsSecion
      elseif occursin("Status", line)
        sec = StatusSection
      elseif occursin("Modelling time", line)
        sec = ModellingTimeSection
      elseif occursin("Computation time", line)
        sec = ComputationTimeSection
      elseif occursin("STAN version", line)
        sec = VersionSection
      elseif occursin("Machine information", line)
        sec = MachineInfosSection
      else
        if sec == NEqsSection
          neqs = parse(Int64, line)
          sec = nothing
        elseif sec == NVarsSecion
          nvars = parse(Int64, line)
          sec = nothing
        elseif sec == StatusSection
          stan_status = line
          sec = nothing
        elseif sec == ModellingTimeSection
          modelling_time = parse(Float64, line)
          sec = nothing
        elseif sec == ComputationTimeSection
          computation_time = parse(Float64, line)
          sec = nothing
        elseif sec == VersionSection
          version = line
          sec = nothing
        elseif sec == MachineInfosSection
          key, val = split(line, ":")
          set!(machine_infos, key, String(strip(val)))
        end
      end
    end
  end
  if iserror(fileresult)
    throw(FileIOError(fileresult.msg))
  end

  STANInfos(neqs, nvars, stan_status, modelling_time, computation_time, version, machine_infos)
end

function add_STAN_infos!(infos::Dictionary{String, Any}, filepath::String)
  try
    stan_infos = parse(STANInfos, filepath)

    add_info!(infos, INFO_STAN_STATUS, stan_infos.stan_status)
    add_info!(infos, INFO_STAN_MODELLING_TIME, stan_infos.modelling_time)
    add_info!(infos, INFO_STAN_COMPUTATION_TIME, stan_infos.computation_time)
    add_info!(infos, INFO_STAN_VERSION, stan_infos.version)
    add_info!(infos, INFO_STAN_MACHINE, stan_infos.machine_infos)
  catch
    @warn "Reading STAN solving details wasn't successful."
  end
end

# ====================
# STAN Trace parsing
# ====================

"""
Struct which holds result data from STAN
"""
struct STANResult
  "Values assigned to variables"
  variables_value::Dictionary{Label, Float64}
  "Unassigned variables"
  unassigned_variables::Dictionary{Label, UnmeasuredVariable}
  "Uncertainties assigned to variables"
  uncertainties::Dictionary{Label, Float64}
end

function STANResult(;
  variables_value::Dictionary{Label, Float64} = Dictionary{Label, Float64}(),
  unassigned_variables::Dictionary{Label, UnmeasuredVariable} = Dictionary{Label, UnmeasuredVariable}(),
  uncertainties::Dictionary{Label, Float64} = Dictionary{Label, Float64}(),
)
  STANResult(variables_value, unassigned_variables, uncertainties)
end

"""
Struct which holds data from reconciliation optimization produced by STAN.
"""
@kwdef struct STANTrace
  "Model's name."
  name::String = ""
  "Variables at the initial (pre-reconciliation) state."
  initial_variables::Vector{Variable}
  "Reconciliation results at the initial state."
  initial_results::STANResult
  "If the model is infeasible we might get no iterations in the trace."
  number_iterations::Union{Nothing, Int64}
  "Mapping from equation identifiers to their string descriptions."
  equations_legend::Dictionary{String, String}
  "Equations of the reconciliation problem."
  equations::Vector{Equation}
  "Reconciliation results at the final (converged) state."
  final_results::STANResult
  "Reconciliation results using the Kelly method, or `missing` if not available."
  kelly_results::Union{STANResult, Missing}
  infos::Dictionary{String, Any} = Dictionary{String, Any}()
end

function Base.show(io::IO, trace::STANTrace)
  nvariables = length(trace.initial_variables)
  nequations = length(trace.equations)
  println(io, "• STAN Trace with:")
  println(io, "  • Variables: $nvariables")
  println(io, "  • Equations: $nequations")
end

"""
Remove trailing comma from a string if present.
"""
strip_trailing_comma(s::AbstractString)::String = last(s) == ',' ? s[1:(end - 1)] : s

"""
Helper struct to hold parsed variable line data before storage.
"""
struct ParsedVariableLine
  "Variable name (key)."
  name::String
  "Additional information or description associated with the variable."
  info::String
  "Whether the variable represents a transfer coefficient."
  tc::Bool
  "Kind of variable value: `:unmeasured`, `:measured`, or `:fixed`."
  value_type::Symbol
  "Parsed numeric value, or `nothing` for unmeasured variables."
  value::Union{Float64, Nothing}
  "Parsed uncertainty, or `nothing` if not available."
  uncertainty::Union{Float64, Nothing}
end

"""
Parse a variable line in STAN format: "Value | Key | Info"

Returns a `ParsedVariableLine` containing the parsed components.
The parsing handles three formats:
- "?" for unmeasured variables
- "value±uncertainty" for measured variables
- "value" for fixed variables
"""
function parse_variable_line(line::String)::ParsedVariableLine
  val, key, info = split(line, "|")

  # Remove spaces and trailing commas.
  value_str = replace(strip(val), "," => ".")
  name = String(strip(key))
  info = strip_trailing_comma(String(strip(info)))

  # Transfer coefficients are parsed from `info`.
  tc = occursin("Transfer coefficient", info) || occursin("Transferkoeffizient", info)

  # Parse value based on format.
  if value_str == "?"
    ParsedVariableLine(name, info, tc, :unmeasured, nothing, nothing)
  elseif occursin("±", value_str)
    value, uncertainty = split(value_str, "±")
    value = parse(Float64, value)
    uncertainty = parse(Float64, strip(uncertainty))
    ParsedVariableLine(name, info, tc, :measured, value, uncertainty)
  else
    value = parse(Float64, value_str)
    ParsedVariableLine(name, info, tc, :fixed, value, nothing)
  end
end

function parse_variable!(target::Vector{Variable}, line::String)::Vector{Variable}
  parsed = parse_variable_line(line)

  var = if parsed.value_type == :unmeasured
    UnmeasuredVariable(; name = parsed.name, info = parsed.info, tc = parsed.tc)
  elseif parsed.value_type == :measured
    value = parsed.value
    uncertainty = parsed.uncertainty
    name = parsed.name
    info = parsed.info
    tc = parsed.tc

    if iszero(uncertainty)
      # Exception: treat entries like
      # 10±0            | FV:292               | F25:Mass flow:N:2013
      # which have zero uncertainy as fixed, since otherwise the weights (Q) matrix is non-invertible.
      FixedVariable(; name, value, info, tc)
    else
      if !relaxed_geq(value, 0.0)
        @warn "Mean of measured variable $name is smaller than 0.0 - dropping lower bound."
        lb = -Inf
      else
        lb = missing
      end
      if tc && !relaxed_leq(value, 1.0)
        @warn "Mean of measured transfer coefficient $name is smaller than 1.0 - dropping upper bound."
        ub = Inf
      else
        ub = missing
      end
      MeasuredVariable(; name, value, uncertainty, info, tc, lb, ub)
    end
  else # :fixed
    FixedVariable(; name = parsed.name, value = parsed.value, info = parsed.info, tc = parsed.tc)
  end
  push!(target, var)
end

function parse_variable_value!(target::STANResult, line::String)::STANResult
  parsed = parse_variable_line(line)

  if parsed.value_type == :unmeasured
    set!(
      target.unassigned_variables,
      parsed.name,
      UnmeasuredVariable(; name = parsed.name, info = parsed.info, tc = parsed.tc),
    )
  elseif parsed.value_type == :measured
    set!(target.variables_value, parsed.name, parsed.value)
    set!(target.uncertainties, parsed.name, parsed.uncertainty)
  else # :fixed
    set!(target.variables_value, parsed.name, parsed.value)
  end
  target
end

function add_equation_legend!(target::Dictionary{String, String}, line::String)::Dictionary{String, String}
  key, info = split(line, "|")
  # Remove spaces and trailing commas.
  key = strip(key)
  info = strip_trailing_comma(strip(info))
  set!(target, key, info)
end

function parse_iterations(line::String)::Int64
  _, cnt = split(line, ":")
  parse(Int64, cnt)
end

"Helper to parse an equation from STAN given as a string."
function parse_equation!(equations, line)
  parsed_eq = parse_equation(replace(line, "," => "."))
  push!(equations, parsed_eq)
end

"Parse Kelly's STAN trace."
function parse_kelly_trace(path::String)
  # Initialize the dictionary to store the trace.
  kelly_results = STANResult()

  open(path, "r") do file
    for line in Iterators.drop(eachline(file), 2)
      line = String(strip(line))

      if line == "[" || line == "]" || isempty(line)
        continue
      end

      parse_variable_value!(kelly_results, line)
    end

  end
  return kelly_results
end

@enum Section InitialVariablesSection InitialResultsSection LegendSection EquationsSection IterationsSection FinalResultsSection

"""
Section headers in English.
"""
const section_header = dictionary([
  InitialVariablesSection => "Vector of variables (Value | Key | Info) = ",
  InitialResultsSection => "Vector of results (Value | Key | Info) = ",
  LegendSection => "Legend (Key | Info):",
  EquationsSection => "Equations:",
  IterationsSection => "Number of iterations:",
  FinalResultsSection => "Vector of results (Value | Key | Info) = ",
])

"""
Section headers in German.
"""
const section_header_DE = dictionary([
  InitialVariablesSection => "Variablenvektor (Value | Key | Info) = ",
  InitialResultsSection => "Ergebnisvektor (Value | Key | Info) = ",
  LegendSection => "Legende (Key | Info):",
  EquationsSection => "Gleichungen:",
  IterationsSection => "Anzahl der Iterationen:",
  FinalResultsSection => "Ergebnisvektor (Value | Key | Info) = ",
])

"""
Check if a line matches a section header in either English or German.
"""
function matches_section_header(line::String, section::Section)::Bool
  line == section_header[section] || line == section_header_DE[section]
end

"""
Check if a line contains a section header in either English or German.
"""
function contains_section_header(line::String, section::Section)::Bool
  occursin(section_header[section], line) || occursin(section_header_DE[section], line)
end

# TODO Define safe_parse and use it here.
function Base.parse(::Type{STANTrace}, path::String)::STANTrace
  filename = joinpath(path, "data", "trace.txt")
  name = basename(path)
  if !isfile(filename)
    throw(ArgumentError("Expected trace file under `/data/trace.txt`."))
  end

  initial_variables = Variable[]
  initial_results = STANResult()
  number_iterations = nothing
  equations_legend = Dictionary{String, String}()
  equations = Equation[]
  final_results = STANResult()

  fileresult = open(filename, "r") do file
    previous = nothing
    current = nothing

    for line in eachline(file)
      if isnothing(current)
        if isnothing(previous) && matches_section_header(line, InitialVariablesSection)
          current = InitialVariablesSection
        elseif previous == InitialVariablesSection && matches_section_header(line, InitialResultsSection)
          current = InitialResultsSection
        elseif previous == InitialResultsSection && matches_section_header(line, LegendSection)
          current = LegendSection
        elseif previous == InitialResultsSection && matches_section_header(line, EquationsSection)
          # Some files do not have a legends section and go directly to equations.
          current = EquationsSection
        elseif previous == LegendSection && matches_section_header(line, EquationsSection)
          current = EquationsSection
        elseif previous == EquationsSection && contains_section_header(line, IterationsSection)
          number_iterations = parse_iterations(line)
          previous = IterationsSection
          current = nothing
          continue
        elseif matches_section_header(line, FinalResultsSection) &&
               # In some cases (eg. Status Quo) the iterations section is missing.
               # But the presence of equations section indicates that we are now parsing the final results vector.
               (previous == IterationsSection || previous == EquationsSection)
          current = FinalResultsSection
        end
        !isnothing(current) && continue
      end

      if line == "["
        continue
      elseif line == "]"
        previous = current
        current = nothing
        continue
      end
      if current == InitialVariablesSection
        parse_variable!(initial_variables, line)
      elseif current == InitialResultsSection
        parse_variable_value!(initial_results, line)
      elseif current == LegendSection
        add_equation_legend!(equations_legend, line)
      elseif current == EquationsSection
        if isempty(line)
          previous = current
          current = nothing
          continue
        end
        parse_equation!(equations, line)
      elseif current == FinalResultsSection
        parse_variable_value!(final_results, line)
      end
    end
  end
  if iserror(fileresult)
    throw(FileIOError(fileresult.msg))
  end

  if isempty(equations)
    throw(
      ArgumentError(
        "No equations found, make sure that `Display Equations` from STAN's Calculation Dialog is activated before saving the trace.",
      ),
    )
  end

  kelly_file = joinpath(path, "data", "trace_Kelly.txt")
  kelly_results = isfile(kelly_file) ? parse_kelly_trace(kelly_file) : missing

  infos = Dictionary{String, Any}()
  add_STAN_infos!(infos, path)

  STANTrace(;
    name,
    initial_variables,
    initial_results,
    number_iterations,
    equations_legend,
    equations,
    final_results,
    kelly_results,
    infos,
  )
end

end
