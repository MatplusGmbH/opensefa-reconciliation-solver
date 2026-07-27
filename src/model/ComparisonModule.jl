"""
Convenience functions to compare numerical values up to some user-controlled precision.

Reference implementation: https://github.com/JuliaReach/ReachabilityBase.jl/tree/master
"""
module ComparisonModule

using Reexport
using Dictionaries

export
  set_rtol,
  set_ztol,
  set_atol,
  set_tolerance,
  get_rtol,
  get_ztol,
  get_atol,
  default_tolerance,
  relaxed_geq,
  relaxed_leq,
  relaxed_isapprox,
  isapprox_dict,
  isapproxzero,
  tight_geq,
  tight_leq

"""
Type that represents the tolerances for a given numeric type.

### Notes

The type `Tolerance`, parametric in the numeric type `N`, is used to store default
values for numeric comparisons. It is mutable and setting the value of a field
affects the getter functions hence it can be used to globally fix the tolerance.

Default values are defined for the most commonly used numeric types, and for those
cases when other numeric types are needed one can extend the default values
as explained next.

The cases `Float64` and `Rational` are special in the sense that they are the most
commonly used types in applications. Getting and setting default tolerances
is achieved with the functions `get_rtol` and `set_rtol` (and similarly for the other
tolerances); the implementation creates an instance of `Tolerance{Float64}`
(resp. `Tolerance{Rational}`) and sets some default values. Since `Tolerance`
is mutable, setting a value is possible e.g. `set_rtol(Type{Float64}, ε)` for some
floating-point `ε`.

For all other cases, a dictionary mapping numeric types to instances of `Tolerance`
for that numeric type is used. For floating-point types, a default value has been
defined through `default_tolerance` as follows:

```julia
default_tolerance(N::Type{<:AbstractFloat}) = Tolerance(Base.rtoldefault(N), N(10) * sqrt(eps(N)), zero(N))
```
Hence to set a single tolerance (either `rtol`, `ztol` or `atol`) for a given
floating-point type, use the corresponding `set_rtol` function, while the values
which have not been set will be pulled from `default_tolerance`. If you would like
to define the three default values at once, or are computing with a non floating-point
numeric type, you can just extend `default_tolerance(N::Type{<:Number})`.
"""
mutable struct Tolerance{N <: Number}
  "Relative tolerance."
  rtol::N
  "Zero tolerance or absolute tolerance for comparison against zero."
  ztol::N
  "Absolute tolerance."
  atol::N
end
default_tolerance(N::Type{<:Rational}) = Tolerance(zero(N), zero(N), zero(N))
default_tolerance(N::Type{<:Integer}) = Tolerance(zero(N), zero(N), zero(N))
default_tolerance(N::Type{<:AbstractFloat}) = Tolerance(Base.rtoldefault(N) * N(10_000), N(10) * sqrt(eps(N)), zero(N))

function set_tolerance(N, tolerance::Tolerance = default_tolerance(N))
  set_rtol(N, tolerance.rtol)
  set_ztol(N, tolerance.ztol)
  set_atol(N, tolerance.atol)
  return tolerance
end

# Global Float64 tolerances.
const DEFAULT_TOL_FLOAT64 = default_tolerance(Float64)

get_rtol(::Type{Float64}) = DEFAULT_TOL_FLOAT64.rtol
get_ztol(::Type{Float64}) = DEFAULT_TOL_FLOAT64.ztol
get_atol(::Type{Float64}) = DEFAULT_TOL_FLOAT64.atol

set_rtol(::Type{Float64}, ε::Float64) = DEFAULT_TOL_FLOAT64.rtol = ε
set_ztol(::Type{Float64}, ε::Float64) = DEFAULT_TOL_FLOAT64.ztol = ε
set_atol(::Type{Float64}, ε::Float64) = DEFAULT_TOL_FLOAT64.atol = ε

# Global rational tolerances.
const DEFAULT_TOL_RATIONAL = default_tolerance(Rational)

get_rtol(::Type{<:Rational}) = DEFAULT_TOL_RATIONAL.rtol
get_ztol(::Type{<:Rational}) = DEFAULT_TOL_RATIONAL.ztol
get_atol(::Type{<:Rational}) = DEFAULT_TOL_RATIONAL.atol

set_rtol(::Type{<:Rational}, ε::Rational) = DEFAULT_TOL_RATIONAL.rtol = ε
set_ztol(::Type{<:Rational}, ε::Rational) = DEFAULT_TOL_RATIONAL.ztol = ε
set_atol(::Type{<:Rational}, ε::Rational) = DEFAULT_TOL_RATIONAL.atol = ε

# Global default tolerances for other numeric types.
const DEFAULT_TOL_N = Dict{Type{<:Number}, Tolerance}()

get_rtol(N::Type{<:Number})::N = get!(DEFAULT_TOL_N, N, default_tolerance(N)).rtol
get_ztol(N::Type{<:Number})::N = get!(DEFAULT_TOL_N, N, default_tolerance(N)).ztol
get_atol(N::Type{<:Number})::N = get!(DEFAULT_TOL_N, N, default_tolerance(N)).atol

function set_rtol(N::Type{NT}, ε::NT) where {NT <: Number}
  tol = get!(DEFAULT_TOL_N, N, default_tolerance(N))
  tol.rtol = ε
  ε
end

function set_ztol(N::Type{NT}, ε::NT) where {NT <: Number}
  tol = get!(DEFAULT_TOL_N, N, default_tolerance(N))
  tol.ztol = ε
  ε
end

function set_atol(N::Type{NT}, ε::NT) where {NT <: Number}
  tol = get!(DEFAULT_TOL_N, N, default_tolerance(N))
  tol.atol = ε
  ε
end

"Define isapprox for Dictionaries."
function isapprox_dict(dict1::Dictionary, dict2::Dictionary; kwargs...)::Bool
  if !issetequal(keys(dict1), keys(dict2))
    return false
  end
  perm = indexin(collect(keys(dict2)), collect(keys(dict1)))
  isapprox(collect(dict1)[perm], collect(dict2); kwargs...)
end

"""
Determine if `x` is approximately zero.

It is considered that `x ≈ 0` whenever `x` (in absolute value) is smaller
than the tolerance for zero, `ztol`.
"""
function isapproxzero(x::N; ztol::Real = get_ztol(N)) where {N <: Real}
  abs(x) <= ztol
end

"""
Determine if `x` is approximately equal to `y`.

### Algorithm

We first check if `x` and `y` are both approximately zero, using
`isapproxzero(x) && isapproxzero(y)`.
If that fails, we check if `x ≈ y`, using Julia's `isapprox(x, y)`.
In the latter check we use `atol` absolute tolerance and `rtol` relative
tolerance.

Comparing to zero with default tolerances is a special case in Julia's
`isapprox`, see the last paragraph in `?isapprox`. This function tries to
combine `isapprox` with its default values and a branch for `x ≈ 0 ≈ y` which
includes `x == 0 == y` but also admits a tolerance `ztol`.

Note that if `x = ztol` and `y = -ztol`, then `|x-y| = 2*ztol` and still
`relaxed_isapprox` returns `true`.
"""
function relaxed_isapprox(x::N, y::N;
  rtol::Real = get_rtol(N),
  ztol::Real = get_ztol(N),
  atol::Real = get_atol(N)) where {N <: Real}
  if isapproxzero(x; ztol = ztol) && isapproxzero(y; ztol = ztol)
    true
  else
    isapprox(x, y; rtol = rtol, atol = atol)
  end
end

# Different numeric types with promotion.
function relaxed_isapprox(x::N, y::M; kwargs...) where {N <: Real, M <: Real}
  relaxed_isapprox(promote(x, y)...; kwargs...)
end

"""
Determine if `x` is smaller than or approximately equal to `y`.

The `x <= y` comparison is split into `x < y` or `x ≈ y`; the latter is
implemented by extending Julia's built-in `isapprox(x, y)` with an absolute
tolerance that is used to compare against zero.
"""
function relaxed_leq(
  x::N,
  y::N;
  rtol::Real = get_rtol(N),
  ztol::Real = get_ztol(N),
  atol::Real = get_atol(N),
) where {N <: AbstractFloat}
  x <= y || relaxed_isapprox(x, y; rtol = rtol, ztol = ztol, atol = atol)
end

# Fallback implementation performs exact comparison.
function relaxed_leq(x::N, y::N; kwargs...) where {N <: Real}
  x <= y
end

# Implementation with type promotion.
function relaxed_leq(x::N, y::M; kwargs...) where {N <: Real, M <: Real}
  relaxed_leq(promote(x, y)...; kwargs...)
end

"Determine if `x` is greater than or equal to `y`."
function relaxed_geq(x::Real, y::Real; kwargs...)
  relaxed_leq(y, x; kwargs...)
end

"Determine if `x` is less than but not approximately equal to `y`."
function tight_leq(x::N, y::N; kwargs...) where {N}
  relaxed_leq(x, y; kwargs...) && !relaxed_isapprox(x, y; kwargs...)
end

"Determine if `x` is greater than but not approximately equal to `y`."
function tight_geq(x::N, y::N; kwargs...) where {N}
  relaxed_geq(x, y; kwargs...) && !relaxed_isapprox(x, y; kwargs...)
end

end
