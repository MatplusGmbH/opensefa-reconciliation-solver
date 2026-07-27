# [Getting Started](@id getting_started)

This guide will help you get up and running with OpenSEFA.jl.

## Installation

To install this package, clone the repository and instantiate it using the Julia package manager:

```julia
$ git clone https://github.com/MatplusGmbH/opensefa-reconciliation-solver.git
$ cd opensefa-reconciliation-solver
$ julia --project
```

In the Julia REPL, activate the package manager mode by typing `]`, then instantiate:

```julia
julia> ]
pkg> instantiate
```

This downloads and precompiles all dependencies. The first time this runs, it may take a few minutes.

!!! note
    The `instantiate` command is only required once, or when dependencies change.

## Basic Workflow

Solving a data reconciliation problem with OpenSEFA.jl involves three main steps:

### 1. Specification

Load your model data into the internal representation (`ReconciliationProblem`). You can do this by:

- Defining problems programmatically in Julia
- Loading from STAN trace files
- Loading from JSON format
- Using the HTTP API

### 2. Solution

Solve the reconciliation problem using default or custom algorithm choices to obtain a `ReconciliationSolution`.

### 3. Validation

Verify the quality of the solution by examining:
- Residual terms magnitude
- Constraint satisfaction
- Objective value

You can iterate between steps 2 and 3, adjusting solver settings as needed.

## Quick Example

Here's a complete example that loads and solves a simple problem:

```julia
# Load the package
using OpenSEFA

# Load the seven flows example from the test suite
include(joinpath(pkgdir(OpenSEFA), "test", "examples", "seven_flows.jl"))

# Parse the equations and load into Julia data structures
rec = seven_flows_cencic_2012()

# Solve the problem with the default algorithm
sol = solve(rec)
```

The `rec` object shows the problem structure:

```julia
julia> rec
• Data reconciliation optimization problem with:
  • Name: Seven flows Cencic2016
  • Equations: 4
  • Variables: 8
    • Measured: 4
    • Unmeasured: 3
    • Fixed: 1
```

The `sol` object contains the solution:

```julia
julia> sol
• Solution of data reconciliation optimization problem with:
  • Name: Seven flows Cencic2016
  • Equations: 4
  • Variables: 8
    • Measured: 4
    • Unmeasured: 3
    • Fixed: 1
  • Objective value: 0.2959
  • Equations discrepancy (max abs): 1.026e-9
  • Fixed variables discrepancy (max abs): 0.0
```

## Loading STAN Models

OpenSEFA.jl can import models from STAN software. There are two main approaches:

### Loading STAN Data (CSV files)

Export your data from STAN's data explorer to CSV files, then:

```julia
using OpenSEFA

path = joinpath(MODELS_PATH, "1a_Example")
sdata = load_stan_data(path)
```

This loads all flows and variable specifications into a `STANData` structure.

### Loading STAN Trace Files

To load a complete assembled model from STAN:

1. Run the STAN solver with the verbose option for printing equations enabled
2. Save all results from the trace console into a `trace.txt` file in your model directory
3. Load and solve:

```julia
using OpenSEFA

path = joinpath(MODELS_PATH, "4_PVC Austria (1950-1994)")
strace = parse(STANTrace, path)
rec = convert(ReconciliationProblem, strace)
sol = solve(rec)
```

## Running Tests

To verify your installation, run the test suite:

```julia
# From the Julia REPL with OpenSEFA activated:
julia> ]
pkg> test
```

Or from the terminal:

```bash
$ julia --project test/runtests.jl
```

## Next Steps

- Explore [Examples](@ref examples_page) for more complex use cases
- Learn about the [HTTP API](@ref http_api) for web integration
- Check the API Reference for detailed function documentation
- Read the [References](@ref refs_page) for theoretical background
