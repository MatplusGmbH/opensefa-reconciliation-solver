using OpenSEFA
using Documenter

DocMeta.setdocmeta!(OpenSEFA, :DocTestSetup, :(using OpenSEFA); recursive = true)

makedocs(;
  modules = [OpenSEFA],
  repo = Documenter.Remotes.GitHub("MatplusGmbH", "opensefa-reconciliation-solver"),
  sitename = "OpenSEFA.jl",
  checkdocs = :none,
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://matplusgmbh.github.io/opensefa-reconciliation-solver/stable",
    edit_link = "main",
    assets = String[]),
  pages = [
    "Home" => "index.md",
    "Getting Started" => "getting-started.md",
    "Solving Strategies" => "solving_strategies.md",
    "Examples" => "examples.md",
    "HTTP API" => "http-api.md",
    "API Reference" => [
      "Model" => [
        "ConstraintsModule" => "lib/model/ConstraintsModule.md",
        "ConstraintParserModule" => "lib/model/ConstraintParserModule.md",
        "ComparisonModule" => "lib/model/ComparisonModule.md",
        "ErrorModule" => "lib/model/ErrorModule.md"],
      "Algorithms" => [
        "ReconciliationModule" => "lib/algorithms/ReconciliationModule.md",
        "PresolveModule" => "lib/algorithms/PresolveModule.md",
        "QCQPModule" => "lib/algorithms/QCQPModule.md",
        "RowEchelonModule" => "lib/algorithms/RowEchelonModule.md",
        "SparseRowEchelonModule" => "lib/algorithms/SparseRowEchelonModule.md"],
      "Solvers" => [
        "DefaultSolverModule" => "lib/solvers/DefaultSolverModule.md",
        "JuMPInterfaceModule" => "lib/solvers/JuMPInterfaceModule.md",
        "NonlinearSolveInterfaceModule" => "lib/solvers/NonlinearSolveInterfaceModule.md",
        "IpoptSolverModule" => "lib/solvers/IpoptSolverModule.md"],
      "I/O" => [
        "STANModule" => "lib/io/STANModule.md",
        "SerializationModule" => "lib/io/SerializationModule.md",
        "ExcelModule" => "lib/io/ExcelModule.md",
        "HTTPServerModule" => "lib/io/HTTPServerModule.md"],
    ],
    "Contributing" => "contributing.md",
    "License" => "license.md",
    "References" => "references.md",
  ])

deploydocs(;
  repo = "github.com/MatplusGmbH/opensefa-reconciliation-solver.git",
  push_preview = true)
