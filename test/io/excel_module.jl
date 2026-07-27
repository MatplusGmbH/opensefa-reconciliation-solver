using Test
using OpenSEFA.STANModule
using OpenSEFA.ReconciliationModule
using OpenSEFA.ExcelModule
using OpenSEFA.DefaultSolverModule
using XLSX: XLSXFile, readxlsx

@testset "Export to excel file." begin
  if !isdir(MODELS_PATH)
    @warn "Skipping excel export test: MODELS_PATH=$(MODELS_PATH) not found."
  else
    # Instantiate problem.
    path = joinpath(MODELS_PATH, "1_SEFMN_Example")
    strace = parse(STANTrace, path)
    rec = convert(ReconciliationProblem, strace)

    # Solve and export solution to excel file.
    sol = solve(rec)
    filename = tempname() * ".xls"
    export_excel(sol, filename)

    # Open the file.
    content = readxlsx(filename)
    @test content isa XLSXFile
    @test content[1].name == "Sheet1"
    @test content[2].name == "1_SEFMN_Example"
  end
end
