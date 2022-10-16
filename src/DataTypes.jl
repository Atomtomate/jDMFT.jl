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

