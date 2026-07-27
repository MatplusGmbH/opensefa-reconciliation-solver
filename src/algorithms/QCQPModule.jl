"""
Reformulate a nonlinear reconciliation problem as a quadratically constrained quadratic
program (QCQP).

The QCQP form follows the notation of Boyd & Vandenberghe (2007):

    minimize  f0(x) = xᵀ P0 x + q0ᵀ x + r0
    s.t.      fi(x) = xᵀ Pi x + qiᵀ x + ri = 0,  i = 1, …, m

where `P0 = Q⁻¹` is the inverse of the measurement covariance matrix. Variables are
ordered as measured → unmeasured → fixed.

See also: https://stanford.edu/~boyd/papers/pdf/qcqp.pdf
"""
module QCQPModule

using Dictionaries
using SparseArrays
using LinearAlgebra
using Reexport

using ..ConstraintsModule
using ..ReconciliationModule

@reexport import ..ConstraintsModule: is_linear

export
  QuadraticConstraint,
  QuadraticallyConstrainedQuadraticProgram,
  QCQP,
  evaluate_unsafe,
  evaluate_gradient_unsafe!

"""
Constraint of the form `x^T * P * x + q^T * x + r = 0`.

If there are no bilinear terms, then `iszero(P)` holds.
"""
struct QuadraticConstraint{VT <: AbstractVector, MT}
  "Quadratic coefficient matrix (symmetric by default)."
  P::MT
  "Linear coefficient vector."
  q::VT
  "Constant term."
  r::Float64
end

"""
Convert a bilinear equation into a quadratic constraint.
"""
function Base.convert(
  ::Type{QuadraticConstraint},
  eq::Equation,
  variable_index::Dictionary{Label, Int64};
  symmetric::Bool = true,
)
  # Evaluate constant value.
  r = eq.constant_term.value

  # Evaluate linear vector.
  nvars = length(variable_index)
  q = spzeros(Float64, nvars)
  for t in eq.linear_terms
    # Accumulation is not needed if the equation has all distinct linear terms
    # (same for the iteration over the bilinear terms above).
    idx = variable_index[t.name]
    q[idx] += t.factor
  end

  # Evaluate bilinear matrix.
  P = spzeros(Float64, nvars, nvars)
  for t in eq.bilinear_terms
    idx1 = variable_index[t.name1]
    idx2 = variable_index[t.name2]

    if symmetric
      # Ensure that P is symmetric.
      P[idx1, idx2] += t.factor / 2
      P[idx2, idx1] += t.factor / 2

    else
      P[idx1, idx2] += t.factor
    end
  end
  P = symmetric ? Symmetric(P) : P
  QuadraticConstraint(P, q, r)
end

"Whether the constraint is linear."
is_linear(cons::QuadraticConstraint) = iszero(cons.P)

"""
Representation of a nonlinear reconciliation optimization problem as a
quadratically constrained quadratic program.
"""
struct QuadraticallyConstrainedQuadraticProgram
  "Original problem."
  problem::ReconciliationProblem
  "Map variables to the corresponding index in the system."
  variable_index::Dictionary{Label, Int64}
  "f0(x) = x^T * P0 * x + q0^T * x + r0"
  objective::Union{QuadraticConstraint, Nothing}
  "fi(x) = x^T * Pi * x + qi^T * x + ri = 0 for all i = 1, ..., m"
  constraints::Vector{QuadraticConstraint}
  "Vector of lower bounds for each variable (missing if it not specified)."
  lower_bounds::Dictionary{Label, Union{Float64, Missing}}
  "Vector of upper bounds for each variable (missing if it not specified)."
  upper_bounds::Dictionary{Label, Union{Float64, Missing}}
end

"Alias for quadratically constrained quadratic program."
const QCQP = QuadraticallyConstrainedQuadraticProgram

"Convert a reconciliation problem into a quadratically constrained quadratic program."
function Base.convert(::Type{QuadraticallyConstrainedQuadraticProgram}, rec::ReconciliationProblem)
  variable_index = Dictionary{Label, Int64}()

  L = length(variable_index)
  for (i, v) in enumerate(rec.measured)
    set!(variable_index, v.name, L + i)
  end

  L = length(variable_index)
  for (i, v) in enumerate(rec.unmeasured)
    set!(variable_index, v.name, L + i)
  end

  L = length(variable_index)
  for (i, v) in enumerate(rec.fixed)
    set!(variable_index, v.name, L + i)
  end

  # Objective function as a quadratic constraint.
  Q = weights_matrix(rec)
  objective = if isnothing(Q)
    nothing
  else
    P0 = inv(Q)
    xm = [x.value for x in rec.measured]
    q0 = -(P0 + P0') * xm
    r0 = dot(xm, P0 * xm)
    QuadraticConstraint(P0, q0, r0)
  end

  constraints = QuadraticConstraint[]
  for eq in rec.equations
    ceq = convert(QuadraticConstraint, eq, variable_index)
    push!(constraints, ceq)
  end

  lower_bounds = Dictionary{Label, Union{Float64, Missing}}()
  upper_bounds = Dictionary{Label, Union{Float64, Missing}}()
  for (n, i) in pairs(variable_index)
    var = get_variable(rec, n)
    lb = lower_bound(var)
    ub = upper_bound(var)
    set!(lower_bounds, n, lb)
    set!(upper_bounds, n, ub)
  end

  QuadraticallyConstrainedQuadraticProgram(rec, variable_index, objective, constraints, lower_bounds, upper_bounds)
end

# ====================================
# Evaluation and residual computation

"Evaluate the constraint `f(x) =  r + <q, x> + <x, Px>` on `x = u`."
function evaluate_unsafe(cons::QuadraticConstraint, u::AbstractVector)
  # This computation is valid both for the linear and for the nonlinear case.
  cons.r + dot(cons.q, u) + dot(u, cons.P, u)
end

"Evaluate the gradient `f(x) =  r + <q, x> + <x, Px>` on `x = u` which is `∇f(x) = q + (P + P')x`."
function evaluate_gradient_unsafe!(grad, cons::QuadraticConstraint, u::AbstractVector)::Vector{Float64}
  copy!(grad, cons.q)
  # grad += P * u
  mul!(grad, cons.P, u, 1.0, 1.0)
  # grad += P' * u
  mul!(grad, transpose(cons.P), u, 1.0, 1.0)
end
function evaluate_gradient_unsafe!(
  grad,
  cons::QuadraticConstraint{VT, MT},
  u::AbstractVector,
)::Vector{Float64} where {VT <: AbstractVector, MT <: Symmetric}
  copy!(grad, cons.q)
  mul!(grad, cons.P, u, 2.0, 1.0)
end

end
