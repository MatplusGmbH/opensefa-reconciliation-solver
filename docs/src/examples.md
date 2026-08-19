# [Examples](@id examples_page)

This page provides practical examples demonstrating the various capabilities of OpenSEFA.jl.

## Basic Example: Seven Flows Problem

The seven flows problem from [[C16]](@ref) is a canonical example in data reconciliation for material flow analysis.

```julia
using OpenSEFA

# Load the problem definition
include(joinpath(MODELS_PATH, "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Solve with default settings
sol = solve(rec)

# Inspect the solution
println("Objective value: ", sol.objective_value)
println("Max constraint violation: ", sol.score.residual)
```

## Using Presolve

The presolve functionality applies symbolic and numeric simplifications to reduce problem size:

```julia
using OpenSEFA

include(joinpath(MODELS_PATH, "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Apply presolve to simplify the problem
pre = unwrap(presolve(rec))

# The presolved problem may have fewer variables/equations
println("Original problem: ", rec)
println("Presolved problem: ", pre.pre)

# Solve the presolved problem
sol = solve(pre.pre)
```

!!! tip
    The default solver (`DefaultSolver`) automatically applies presolve internally. Use explicit presolve when you want to inspect the simplified problem structure.

## Custom Solver Backends

OpenSEFA.jl supports multiple solver backends through the JuMP interface.

### Using the Default Solver

The default solver applies three steps:

1. **Presolve**: Simplify the problem through structural reductions. If fully solved, return immediately.
2. **Starting values**: Compute initial values for unmeasured variables using `NonlinearSolve`. If the problem is fully determined, use these as the final solution.
3. **Optimize**: Solve the full weighted least-squares problem using the JuMP framework (with `Ipopt` as default solver). Starting values from step 2 are used as a warm start when available.

```julia
using OpenSEFA

include(joinpath(MODELS_PATH, "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Explicit use of DefaultSolver (equivalent to solve(rec))
sol = solve(rec, DefaultSolver())
```

### Using Alternative NLP Solvers

You can use any JuMP-compatible solver:

```julia
using OpenSEFA
using MadNLP

include(joinpath(MODELS_PATH, "seven_flows.jl"))
rec = seven_flows_cencic_2012()

# Solve using MadNLP instead of Ipopt
sol = solve(rec, JuMPSolver(optimizer=MadNLP.Optimizer))
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
path = joinpath(MODELS_PATH, "Cencic_2016_STAN")

# Parse the STAN trace file
strace = parse(STANTrace, path)

# Convert to ReconciliationProblem
rec = convert(ReconciliationProblem, strace)

# Solve the problem
sol = solve(rec)

println("Solution for Seven Flows model:")
println("  Objective value: ", sol.objective_value)
println("  Variables: ", length(sol.variables_value))
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
path = joinpath(MODELS_PATH, "Cencic_2016_STAN")

# Load data from STAN CSV exports
sdata = load_stan_data(path)

# Inspect the loaded data
println("Flows data: ", sdata.flows)
println("Variables data: ", sdata.flowvalues)
```

## Working with JSON Models

Models can be defined in JSON format for easy integration with web applications:

```julia
using OpenSEFA

# Load a JSON model
json_path = joinpath(MODELS_PATH, "seven_flows.json")
model = load_model(json_path)
```

## Programmatic Problem Construction

You can build problems programmatically using Julia data structures:

```julia
using OpenSEFA

# Define variables
variables = [
    MeasuredVariable(name = "F1", value = 1.0, uncertainty = 0.1),
    UnmeasuredVariable(name = "F2"),
    MeasuredVariable(name = "F3", value = 2.0, uncertainty = 0.3)
]

# Define equations (mass balance: F1 + F3 - F2 = 0)
equations = [
    Equation(
        ConstantTerm(0.0),
        [LinearTerm("F1", 1.0), LinearTerm("F3", 1.0), LinearTerm("F2", -1.0)],
        BilinearTerm[]
    )
]

# Create the problem and add variables/equations
rec = ReconciliationProblem(name = "Simple Example")
add_variable!(rec, variables)
add_equation!(rec, equations)

# Solve
sol = solve(rec)
```

Equations can alternatively be parsed from a string:
```julia
using OpenSEFA

# Define equation F1 + F3 - F2 = 0
eq = parse_equation("F1 + F3 - F2 = 0")
```

## Inspecting Solutions

The solution object provides various methods to inspect results:

```julia
using OpenSEFA

include(joinpath(MODELS_PATH, "seven_flows.jl"))
rec = seven_flows_cencic_2012()
sol = solve(rec)

# Access solution values
println("Reconciled variable values:")
display(sol.variables_value)

# Check solution quality
println("\nSolution quality:")
println("  Objective: ", sol.objective_value)
display(sol.score)

# Check if solution is acceptable
if sol.score.residual < 1e-6
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

    if sol.score.residual > 1e-4
        @warn "Constraint violations detected"
    end

catch e
    @error "Solver failed" exception=e
    # Fall back to alternative solver or settings
    sol = solve(rec, JuMPSolver(optimizer=AlternativeSolver.Optimizer))
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
using Ipopt

# For large problems, explicit presolve can help
rec = # ... large problem ...
presolve_result = presolve(rec)
pre = unwrap(presolve_result).pre

# Solve with custom tolerances (if using JuMP interface)
sol = solve(pre, JuMPSolver(
    optimizer = Ipopt.Optimizer,
    tol = 1e-6,
    max_iter = 1000
))
```

## Next Steps

- Explore the [HTTP API](@ref http_api) for web-based integration
- Check the API Reference for detailed function documentation
- Review the [References](@ref refs_page) for theoretical background
