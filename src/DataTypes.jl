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

