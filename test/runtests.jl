using Test
using SafeTestsets: @safetestset

TEST_FILES = [
  # model/
  "model/error_module.jl",
  "model/comparison_module.jl",
  "model/constraints_module.jl",
  "model/constraint_parser_module.jl",
  # algorithms/
  "algorithms/row_echelon_module.jl",
  "algorithms/sparse_row_echelon_module.jl",
  "algorithms/qcqp_module.jl",
  "algorithms/reconciliation_module.jl",
  "algorithms/presolve_module.jl",
  # solvers/
  "solvers/jump_interface_module.jl",
  "solvers/ipopt_solver_module.jl",
  "solvers/nonlinear_solve_interface_module.jl",
  "solvers/default_solver_module.jl",
  # io/
  "io/stan_module.jl",
  "io/serialization_module.jl",
  "io/excel_module.jl",
  "io/http_server_module.jl",
]

function test(files::Vector)
  @testset "OpenSEFA.jl" begin
    foreach(test, files)
  end
end
function test(file::String)
  @info "Testing $file..."
  path = joinpath(@__DIR__, file)
  @eval @time @safetestset $file begin
    include($path)
  end
end

test([TEST_FILES])
