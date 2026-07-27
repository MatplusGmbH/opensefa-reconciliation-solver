# OpenSEFA.jl

```@meta
CurrentModule = OpenSEFA
```

`OpenSEFA.jl` is a Julia package designed for modeling and solving **nonlinear data reconciliation problems** within the context of **Substance & Energy Flow Analysis (SEFA)** and **Material Flow Analysis (MFA)**. These methodologies serve as powerful decision-support tools across diverse fields, including but not limited to:

- **Environmental management**
- **Waste resource management**
- **Policy assessment**

The primary goal of MFA is to **quantify all flows and stocks** in a defined system, along with their associated uncertainties. For a comprehensive introduction to the theory and applications of MFA, refer to the standard handbook [[BR16]](@ref).

This project has been highly inspired by [STAN Software](https://www.stan2web.net/).

## Key Features

- **Nonlinear optimization**: Solve data reconciliation problems with both linear and bilinear constraints
- **Uncertainty quantification**: Integrates measured variable uncertainties directly into the optimization
- **Customizable modeling**: Supports various types of variables and constraints to adapt to a wide range of MFA scenarios
- **Efficient implementation**: Built on Julia for high performance in numerical computation
- **HTTP API**: RESTful API for integration with external tools and web applications
- **STAN compatibility**: Import and solve models from STAN software

## Getting Started

New to OpenSEFA.jl? Check out the [Getting Started](@ref getting_started) guide.

For practical examples, see the [Examples](@ref examples_page) section.

## [Problem Formulation](@id problem_formulation)

This package addresses the following **nonlinear data reconciliation optimization problem**:

```math
\begin{aligned}
\min_{x,y} \quad & (x - x_m)^T Q^{-1} (x - x_m) \\
\text{s.t.} \quad & f_i(x,y,z) := c_i + a_i^T \begin{bmatrix} x \\ y \\ z \end{bmatrix} + \begin{bmatrix} x \\ y \\ z \end{bmatrix}^T B_i \begin{bmatrix} x \\ y \\ z \end{bmatrix} = 0, \quad i = 1, \dots, k\\
& lb \leq \begin{bmatrix} x \\ y \end{bmatrix} \leq ub \quad \text{(componentwise)}
\end{aligned}
```

### Variables and Parameters

- **Measured variables** ($x$): Known variables with specified mean values ($x_m$) and uncertainties
- **Unmeasured variables** ($y$): Unknown variables to be estimated through optimization
- **Fixed variables** ($z$): Known variables without associated uncertainties
- **Weights matrix** ($Q$): A diagonal matrix whose entries reflect the variances of the measured variables
- **Variable bounds** ($lb$, $ub$): Lower and upper bounds that can be set for every measured and unmeasured variable. If a bound isn't provided (or set to `missing`), a lower bound is per default set to $0.0$ for all variables and for transfer coefficients the upper bound is set to $1.0$ per default.
- **Constraints**: Each constraint $f_i$, $i \in \{1, \dots, k\}$, can be linear (where $B_i$ is the zero matrix) or nonlinear (bilinear)

Each constraint is parameterized by:
-  $c_i$: A scalar constant
-  $a_i$: A vector associated with the linear terms
-  $B_i$: A matrix for bilinear terms

## Documentation Contents

```@contents
Pages = [
    "getting-started.md",
    "examples.md",
    "http-api.md",
    "contributing.md",
    "license.md",
    "references.md"
]
Depth = 2
```

## API Reference

For detailed API documentation, see the API Reference section in the sidebar.

## Resources

- **Source Code**: [GitHub Repository](https://github.com/MatplusGmbH/opensefa-reconciliation-solver)
- **Issue Tracker**: [GitHub Issues](https://github.com/MatplusGmbH/opensefa-reconciliation-solver/issues)
- **License**: [PolyForm Noncommercial 1.0.0](@ref license)
- **STAN Software**: [www.stan2web.net](https://www.stan2web.net/)

## References

See the [References](@ref refs_page) section for a complete list of citations.
