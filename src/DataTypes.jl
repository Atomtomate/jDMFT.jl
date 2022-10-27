"""
    SampleMatrix

Holds samples data for f(τ_ij) in each matrix entry. f depends on the algrithm (Weiss GF for CT-INT, Δ for CT_HYB).
This should be used in conjunction with the config struct (wich saves the τ_ij and additional data) for the specific sampler.

Fields
-------------
- **`data`**        : `Matrix{ComplexF64}`, Data for sample matrix
- **`rowCache`**    : `Vector{ComplexF64}`, Row cache for potential new size
- **`colCache`**    : `Vector{ComplexF64}`, Column cache for potential new size
- **`S`**           : `Complex64`, Cache for bottom right entry (remaining new value besides row and colum cache)
- **`N`**           : `Int`, Used matrix dimension, this differs from the actual data size, due to caching, see also [`update_cache!`](@ref)
"""
mutable struct SampleMatrix
    data::Matrix{ComplexF64}
    rowCache::Vector{ComplexF64}
    colCache::Vector{ComplexF64}
    S::ComplexF64
    N::Int
    function SampleMatrix()
        SampleMatrix(10)
    end
    function SampleMatrix(N::Int)
        new(Matrix{ComplexF64}(undef, N,N), Vector{ComplexF64}(undef, N),
            Vector{ComplexF64}(undef, N), 0 + 0im, 0)
    end
end

"""
    fast_update_incr!(M::SampleMatrix, N::Int)

Compute fast inverse update, using Woodbury matrix idenity: 
    `M_inv = inv(M)`
    `Mnew_inv = M_inv` after inserting row/col/corner obtained from `M.rowCach`, `M.colCache`, `M.S`
    `Mnew = inv(Mnew_inv)`

`M.rowCache`, `M.colCache` and `M.S` must be set to the new values, before calling this function!
`N` is the new size. TODO: currently only rank 1 updates supported, so `N = M.N+1`.
"""
function fast_update_incr!(M::SampleMatrix, N::Int)
    Nold = M.N
    Mi   = view(M.data, 1:Nold, 1:Nold)
    Ri   = transpose(view(M.rowCache, 1:Nold))
    Ci   = view(M.colCache, 1:Nold)
    S    = Nold != 0 ? 1/(M.S - Ri*(Mi * Ci)) : 1/M.S

    #TODO: this can probably be done very fast in a nested loop
    M.data[1:Nold,N]  = -(Mi * Ci) * S
    M.data[N,1:Nold]  = -S .* (Ri * Mi)
    M.data[1:Nold,1:Nold]   = Mi + (Mi*Ci)*S*(Ri*Mi)
    M.data[N,N] =  S
    M.N = N
end

"""
    fast_update_decr!(M::SampleMatrix, N::Int)

Compute fast inverse update, using Woodbury matrix idenity: 
    `M_inv = inv(M)`
    `Mnew_inv = ``M_inv` after deleting right column and bottom row
    `Mnew = inv(Mnew_inv)`

`N` is the new size. TODO: currently only rank 1 updates supported, so `N = M.N-1`.
"""
function fast_update_decr!(M::SampleMatrix, N::Int)
    Ri = view(M.data, M.N, 1:N)
    Ci = view(M.data, 1:N, M.N)
    #TODO: this can probably be done very fast in a nested loop
    M.data[1:N,1:N] = M.data[1:N,1:N] .- (Ci * Ri') / M.data[M.N,M.N]
    M.N = N
end
