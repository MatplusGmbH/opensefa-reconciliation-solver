"""
Reduced row echelon form computation optimized for sparse matrices.

This module extends the functionality of `RowEchelonModule` by adding specialized methods
for sparse matrices that avoid conversion to dense format.
"""
module SparseRowEchelonModule

using LinearAlgebra
using SparseArrays

export sparse_rref, sparse_rref!, sparse_rref_with_pivots, sparse_rref_with_pivots!
export nullspace_sparse, rank_sparse

"""
Compute the reduced row echelon form of a sparse matrix A in-place.
Optimized to maintain sparsity throughout the computation.

Parameters:
- `A`: Input sparse matrix (modified in-place)
- `ɛ`: Tolerance for determining zero pivots

Returns:
- The modified matrix A in RREF.
"""
function sparse_rref!(A::SparseMatrixCSC{T}, ɛ = eps(real(T)) * size(A, 2)) where {T <: Number}
  nr, nc = size(A)
  i = j = 1

  # Pre-allocate arrays for column operations
  buffer = zeros(T, nr)

  while i <= nr && j <= nc
    # Find max absolute value in column j, starting from row i.
    # Only search through non-zero entries.
    max_val = zero(real(T))
    mi = i

    # Extract non-zero entries in this column segment.
    for k in nzrange(A, j)
      row = A.rowval[k]
      if row >= i && abs(A.nzval[k]) > max_val
        max_val = abs(A.nzval[k])
        mi = row
      end
    end

    if max_val <= ɛ
      # No suitable pivot found, move to next column.
      j += 1
      continue
    end

    # Swap rows if needed.
    if mi != i
      # Row swap for sparse matrices: instead of swapping elements directly, we'll swap rows when needed.
      for col in 1:nc
        A[i, col], A[mi, col] = A[mi, col], A[i, col]
      end
    end

    # Get the pivot value.
    pivot = A[i, j]

    # Normalize the pivot row.
    for col in 1:nc
      if !iszero(A[i, col])
        A[i, col] /= pivot
        # Clean up small values due to numerical precision.
        if abs(A[i, col]) < ɛ
          A[i, col] = zero(T)
        end
      end
    end

    # Zero out the rest of the column using the pivot row.
    for row in 1:nr
      if row != i && !iszero(A[row, j])
        factor = A[row, j]

        # Subtract the pivot row multiplied by the factor.
        for col in 1:nc
          if !iszero(A[i, col])
            A[row, col] -= factor * A[i, col]
            # Clean up small values.
            if abs(A[row, col]) < ɛ
              A[row, col] = zero(T)
            end
          end
        end
      end
    end

    i += 1
    j += 1
  end

  # Remove numerical noise and maintain sparsity.
  dropzeros!(A)

  return A
end

"""
Compute the reduced row echelon form of a sparse matrix A in-place, also returning the pivot columns.

Parameters:
- `A`: Input sparse matrix (modified in-place)
- `ɛ`: Tolerance for determining zero pivots

Returns:
- Tuple containing (modified matrix A in RREF, pivot columns).
"""
function sparse_rref_with_pivots!(A::SparseMatrixCSC{T}, ɛ = eps(real(T)) * size(A, 2)) where {T <: Number}
  nr, nc = size(A)
  i = j = 1
  pivots = Int64[]

  while i <= nr && j <= nc
    # Find max absolute value in column j, starting from row i.
    max_val = zero(real(T))
    mi = i

    for k in nzrange(A, j)
      row = A.rowval[k]
      if row >= i && abs(A.nzval[k]) > max_val
        max_val = abs(A.nzval[k])
        mi = row
      end
    end

    if max_val <= ɛ
      # No suitable pivot found, move to next column.
      j += 1
      continue
    end

    # Swap rows if needed
    if mi != i
      for col in 1:nc
        A[i, col], A[mi, col] = A[mi, col], A[i, col]
      end
    end

    # Get the pivot value.
    pivot = A[i, j]

    # Normalize the pivot row.
    for col in 1:nc
      if !iszero(A[i, col])
        A[i, col] /= pivot
        if abs(A[i, col]) < ɛ
          A[i, col] = zero(T)
        end
      end
    end

    # Zero out the rest of the column using the pivot row.
    for row in 1:nr
      if row != i && !iszero(A[row, j])
        factor = A[row, j]

        for col in 1:nc
          if !iszero(A[i, col])
            A[row, col] -= factor * A[i, col]
            if abs(A[row, col]) < ɛ
              A[row, col] = zero(T)
            end
          end
        end
      end
    end

    push!(pivots, j)
    i += 1
    j += 1
  end

  # Remove numerical noise and maintain sparsity.
  dropzeros!(A)

  return A, pivots
end

"""
Compute the reduced row echelon form of a sparse matrix without modifying the input.

Parameters:
- `A`: Input sparse matrix
- `ɛ`: Tolerance for determining zero pivots

Returns:
- A new sparse matrix in RREF
"""
sparse_rref(A::SparseMatrixCSC{T}, ɛ = eps(real(T)) * size(A, 2)) where {T <: Number} = sparse_rref!(copy(A), ɛ)

"""
Compute the reduced row echelon form of a sparse matrix and return pivot columns without modifying the input.

Parameters:
- `A`: Input sparse matrix
- `ɛ`: Tolerance for determining zero pivots

Returns:
- Tuple containing (new matrix in RREF, pivot columns)
"""
sparse_rref_with_pivots(A::SparseMatrixCSC{T}, ɛ = eps(real(T)) * size(A, 2)) where {T <: Number} =
  sparse_rref_with_pivots!(copy(A), ɛ)

"""
Compute the rank of a sparse matrix using the sparse RREF implementation.

Parameters:
- `A`: Input sparse matrix
- `ɛ`: Tolerance for determining zero pivots

Returns:
- The rank of matrix A
"""
function rank_sparse(A::SparseMatrixCSC{T}, ɛ = eps(real(T)) * size(A, 2)) where {T <: Number}
  _, pivots = sparse_rref_with_pivots(A, ɛ)
  return length(pivots)
end

"""
Compute a basis for the nullspace of a sparse matrix using the sparse RREF implementation.

Parameters:
- `A`: Input sparse matrix
- `ɛ`: Tolerance for determining zero pivots

Returns:
- A sparse matrix whose columns form a basis for the nullspace of A.
"""
function nullspace_sparse(A::SparseMatrixCSC{T}, ɛ = eps(real(T)) * size(A, 2)) where {T <: Number}
  R, pivots = sparse_rref_with_pivots(A, ɛ)
  m, n = size(R)

  if length(pivots) == n
    # No free variables, nullspace is empty.
    return spzeros(T, n, 0)
  end

  # Find the free variables (non-pivot columns).
  free_vars = setdiff(1:n, pivots)

  # Create a matrix for the nullspace basis.
  nullspace_basis = spzeros(T, n, length(free_vars))

  # Create a mapping from pivot columns to their row indices.
  pivot_to_row = Dict{Int, Int}()

  # Find which row contains each pivot.
  for j in pivots
    for i in 1:m
      if !iszero(R[i, j]) && abs(R[i, j] - one(T)) < ɛ
        pivot_to_row[j] = i
        break
      end
    end
  end

  # Construct basis vectors for the nullspace.
  for (idx, free_col) in enumerate(free_vars)
    # Set the free variable to 1
    nullspace_basis[free_col, idx] = one(T)

    # Set the dependent variables based on the RREF.
    for pivot_col in pivots
      row = get(pivot_to_row, pivot_col, 0)
      if row > 0
        nullspace_basis[pivot_col, idx] = -R[row, free_col]
      end
    end
  end

  # Clean up small values.
  dropzeros!(nullspace_basis)

  return nullspace_basis
end

end
