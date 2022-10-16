# ================================================================================
#                   Data structures and related functions                        =
# ================================================================================

const τRangeϵ = 0.001

#
#TODO: use more advanced data structure to capture statistics
"""
    Measurements

Accumulation helper, see Gull et al., 2008.

Fields
-------------
- **`samples`**      : `Vector{ComplexF64}`, List of all samples
- **`τGrid`**        : `Vector{Float64}`, sampling points for τ
- **`NSamples`**     : `Int`, number of samples
- **`totalSign`**    : `Int`, total Monte Carlo sign.
"""
mutable struct Measurements
    samples::Vector{ComplexF64}
    τGrid::Vector{Float64}
    NSamples::Int
    totalSign::Int
end

"""
    Measurements(NBins::Int, β::Float64, type::Symbol)

Generates a Measurements struct for a temperature `β` with `NBins` bins and of grid type `type` which may be:
- 'GaussQuad' : generates a Gauss Radau grid and its corresponding τWeights
- 'Riemann'   : generates an evenly distriubted grid where each point has equal weight
"""
function Measurements(NBins::Int, β::Float64, type::Symbol)
    τGrid, τWeights = if type == :GaussQuad
        τGrid_red, τWeights = gaussradau(0, β-τRangeϵ/NBins, NBins)
    elseif type == :Riemann
        τGrid_red, τWeights = riemann(0, β-τRangeϵ/NBins, NBins)
    else
        throw(DomainError(type, "Type must specify pre-implemented τ-grid type (GaussQuad or Riemann)."))
    end
    return Measurements(zeros(ComplexF64, length(τGrid)), τGrid, 0, 0)
end

"""
    CTInt_Confs

Fields
-------------
- **`β`**         : `Float64`, inverse temperature
- **`U`**         : `Float64`, interaction strenght
- **`GWeiss`**    : `τFunction`, noninteracting impurity Green's function
- **`τList`**     : `Array{Float64}`, τ values on which the Green's function is evaluated
- **`τiList`**    : `Array{Float64}`, sampled τ, i-th entry corresponds to i-th row/column in [`SampleMatrix`](@ref)
- **`siList`**    : `Array{Int}`, external ising field, i-th entry corresponds to i-th row/column in [`SampleMatrix`](@ref)
"""
mutable struct CTInt_Confs
    β::Float64
    U::Float64
    GWeiss::τFunction
    τList::Array{Float64}
    τiList::Array{Float64}
    siList::Array{Int}
    CTInt_Confs(GWeiss::τFunction, τList::Vector{Float64}, U::Float64) = 
        CTInt_Confs(GWeiss,τList,U,10)
        function CTInt_Confs(GWeiss::τFunction, τList::Vector{Float64},U::Float64,N::Int)
            new(GWeiss.β, U, GWeiss, τList, Array{Float64}(undef, N), Array{Int}(undef, N))
    end
end
CTInt_Confs(GWeiss::τFunction, U::Float64) =  CTInt_Confs(GWeiss, GWeiss.τGrid, U)


#TODO: use some sort of module instead...
#TODO: move to DataTypes, remove CTInt_confs from arguments
function update_cache!(C::CTInt_Confs, M::SampleMatrix; threshold = 0, incr = 10)
    if M.N - length(M.rowCache) >= threshold
        lNew = size(M.data,1) + incr
        data_cache = similar(M.data, (size(M.data) .+ incr)...)
        data_cache[axes(M.data)...] .= M.data
        M.data = data_cache
        resize!(M.rowCache,lNew)
        resize!(M.colCache,lNew)
        resize!(C.τiList, lNew)
        resize!(C.siList, lNew)
        return 1
    else
        return 0
    end
end

# ================================================================================
#                          Measurement functions                                 =
# ================================================================================

"""
    accumulate!(rng::AbstractRNG, m::Measurements, sign::Int, confs::CTInt_Confs, M::SampleMatrix; with_τshift=true)

Updates Monte Carlo estimate by measuring impurity Green's function for given configuration
stored in `M`.
"""
function accumulate!(rng::AbstractRNG, m::Measurements, sign::Int, 
                  confs::CTInt_Confs, M::SampleMatrix; with_τshift=true)
    τi_list = Vector{Float64}(undef, length(M.rowCache))
    for i in 1:M.N
        τi_list[i] = confs.τiList[i]
        if with_τshift
            τi_list[i] -= rand(rng, Float64) * confs.β
        end
        M.rowCache[i] = draw_Gτ(confs.GWeiss, τi_list[i])
    end
    M.colCache[1:M.N] = view(M.data, 1:M.N, 1:M.N) * view(M.rowCache, 1:M.N)
    N_bins = length(m.samples)
    for i in 1:M.N
        τi = τi_list[i]
        s, ind = τIndex(τi, m.τGrid, confs.β)
        m.samples[ind] += sign * s * M.colCache[i]
    end
    m.NSamples += 1
    m.totalSign += sign
end

#TODO: use convolution theorem here
function measure_GImp_τ(m::Measurements, GWeiss::τFunction)
    GImp_τ = deepcopy(GWeiss.data)
    for i in 1:length(GImp_τ)
        τi = GWeiss.τGrid[i]
        tmp = 0.0
        for (j,τj_sample) in enumerate(m.τGrid)
            tmp += draw_Gτ(GWeiss, τi - τj_sample) * m.samples[j] ./ m.totalSign
        end
        GImp_τ[i] -= tmp
    end
    return GImp_τ
end

#TODO: frequency measurement
