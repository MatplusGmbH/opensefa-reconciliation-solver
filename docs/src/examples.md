# [Examples](@id examples_page)

This page provides practical examples demonstrating the various capabilities of OpenSEFA.jl.

## Basic Example: Seven Flows Problem

The seven flows problem from [[C16]](@ref) is a canonical example in data reconciliation for material flow analysis.

```julia
using OpenSEFA

# Load the problem definition
include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Solve with default settings
sol = solve(rec)

# Inspect the solution
println("Objective value: ", sol.objective_value)
println("Max constraint violation: ", sol.max_equation_discrepancy)
```

## Using Presolve

The presolve functionality applies symbolic and numeric simplifications to reduce problem size:

```julia
using OpenSEFA

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Apply presolve to simplify the problem
pre = presolve(rec)

# The presolved problem may have fewer variables/equations
println("Original problem: ", rec)
println("Presolved problem: ", pre)

# Solve the presolved problem
sol = solve(pre)
```

!!! tip
    The default solver (`DefaultSolver`) automatically applies presolve internally. Use explicit presolve when you want to inspect the simplified problem structure.

## Custom Solver Backends

OpenSEFA.jl supports multiple solver backends through the JuMP interface.

### Using the Default Solver

The default solver applies a two-stage approach:

1. Use `NonlinearSolve` polyalgorithm to find initial values via presolve
2. Refine with `Ipopt` NLP solver

```julia
using OpenSEFA

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Explicit use of DefaultSolver (equivalent to solve(rec))
sol = solve(rec, DefaultSolver())
```

### Using Alternative NLP Solvers

You can use any JuMP-compatible solver:

```julia
using OpenSEFA
using MadNLP

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Solve using MadNLP instead of Ipopt
sol = solve(rec, JuMPSolver(solver=MadNLP.Optimizer))
```

Other compatible solvers include:
- `Ipopt.Optimizer` (default)
- `MadNLP.Optimizer`
- Any other JuMP-compatible NLP solver

## Loading from STAN Software

OpenSEFA.jl can import models created in STAN software.

### Example: Loading STAN Trace File

```julia
using OpenSEFA

# Path to model directory containing trace.txt
path = joinpath(MODELS_PATH, "4_PVC Austria (1950-1994)")

# Parse the STAN trace file
strace = parse(STANTrace, path)

# Convert to ReconciliationProblem
rec = convert(ReconciliationProblem, strace)

# Solve the problem
sol = solve(rec)

println("Solution for PVC Austria model:")
println("  Objective value: ", sol.objective_value)
println("  Variables: ", length(sol.variables))
```

!!! note
    To generate a `trace.txt` file from STAN:
    1. Open your model in STAN
    2. Run the solver with "verbose equations output" enabled
    3. Save the entire trace console output to `trace.txt` in your model directory

### Example: Loading STAN Data (CSV)

```julia
using OpenSEFA

# Path to model directory containing CSV data files
path = joinpath(MODELS_PATH, "1a_Example")

# Load data from STAN CSV exports
sdata = load_stan_data(path)

# Inspect the loaded data
println("Flows data: ", sdata.flows)
println("Variables data: ", sdata.variables)
```

## Working with JSON Models

Models can be defined in JSON format for easy integration with web applications:

```julia
using OpenSEFA
using JSON3

# Load a JSON model
json_path = joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.json")
json_data = read(json_path, String)
model = JSON3.read(json_data)

# Convert to ReconciliationProblem and solve
# (Implementation depends on your JSON schema)
```

## Programmatic Problem Construction

You can build problems programmatically using Julia data structures:

```julia
using OpenSEFA
using Dictionaries

# Define variables
variables = Dictionary(
    ["F1", "F2", "F3"],
    [
        MeasuredVariable(1.0, 0.1),      # value=1.0, uncertainty=0.1
        UnmeasuredVariable(),
        MeasuredVariable(2.0, 0.3)
    ]
)

# Define equations (mass balance: F1 + F3 = F2)
equations = [
    Equation(
        constant_term = 0.0,
        linear_terms = Dictionary(["F1", "F2", "F3"], [1.0, -1.0, 1.0]),
        bilinear_terms = Dictionary()  # No bilinear terms
    )
]

# Create the problem
rec = ReconciliationProblem(
    name = "Simple Example",
    variables = variables,
    equations = equations
)

# Solve
sol = solve(rec)
```

## Inspecting Solutions

The solution object provides various methods to inspect results:

```julia
using OpenSEFA

include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))
rec = seven_flows_cencic_2012()
sol = solve(rec)

# Access solution values
println("Reconciled variable values:")
for (name, value) in sol.variables
    println("  $name = $value")
end

# Check solution quality
println("\nSolution quality:")
println("  Objective: ", sol.objective_value)
println("  Max equation discrepancy: ", sol.max_equation_discrepancy)
println("  Max fixed var discrepancy: ", sol.max_fixed_discrepancy)

# Check if solution is acceptable
if sol.max_equation_discrepancy < 1e-6
    println("  ✓ Solution satisfies constraints")
else
    println("  ✗ Warning: Large constraint violations")
end
```

## Error Handling and Validation

Handle solver failures gracefully:

```julia
using OpenSEFA

rec = # ... your problem definition ...

try
    sol = solve(rec)

    # Validate solution
    if sol.objective_value > 1000.0
        @warn "Large objective value, solution may be poor"
    end

    if sol.max_equation_discrepancy > 1e-4
        @warn "Constraint violations detected"
    end

catch e
    @error "Solver failed" exception=e
    # Fall back to alternative solver or settings
    sol = solve(rec, JuMPSolver(solver=AlternativeSolver.Optimizer))
end
```

## Performance Tips

For large-scale problems:

1. **Use presolve** to reduce problem size
2. **Choose appropriate solvers** - Ipopt is generally robust, MadNLP may be faster for some problems
3. **Provide good initial guesses** when constructing variables
4. **Monitor solver output** to diagnose convergence issues

```julia
using OpenSEFA

# For large problems, explicit presolve can help
rec = # ... large problem ...
pre = presolve(rec)

# Solve with custom tolerances (if using JuMP interface)
sol = solve(pre, JuMPSolver(
    solver = Ipopt.Optimizer,
    solver_options = Dict(
        "tol" => 1e-6,
        "max_iter" => 1000
    )
))
```

## Next Steps

- Explore the [HTTP API](@ref http_api) for web-based integration
- Check the API Reference for detailed function documentation
- Review the [References](@ref refs_page) for theoretical background
