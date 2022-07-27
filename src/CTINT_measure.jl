# ================================================================================
#                   Data structures and related functions                        =
# ================================================================================

#
#TODO: use more advanced data structure to capture statistics
"""
    Measurements

Fields
-------------
- **`samples`**      : `Vector{ComplexF64}`, List of all samples
- **`τWeights`**     : `Vector{Float64}`, weight corresponding to each point on the `τGrid`
- **`τGrid`**        : `Vector{Float64}`, sampling points for τ
- **`NSamples`**     : `Int`, number of samples
- **`totalSign`**    : `Int`, sign
"""

mutable struct Measurements
    samples::Vector{ComplexF64}
    τWeights::Vector{Float64}
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
        τGrid_red, τWeights = gaussradau(NBins);
        τGrid_red .* β ./ 2 .+ β ./ 2, τWeights
    elseif type == :Riemann
        collect(LinRange(0,β-0.001/NBins, NBins)), (ones(Float64, length(NBins)) ./ length(NBins))
    else
        throw(DomainError(type, "Type must specify pre-implemented τ-grid type (GaussQuad or Riemann)."))
    end
    return Measurements(zeros(ComplexF64, length(τGrid)), τWeights, τGrid, 0, 0)
end

"""
    CTInt_Confs

Fields
-------------
- **`β`**         : `Float64`, inverse temperature in units of 1/t
- **`U`**         : `Float64`, interaction strenght in units of t
- **`GWeiss`**    : `τFunction`, noninteracting impurity Green's function
- **`τList`**     : `Array{Float64}`, τ values on which the Green's function is evaluated
- **`τiList`**    : `Array{Float64}`, ???
- **`siList`**    : `Array{Int}`, ???
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
measure_τ!(rng::AbstractRNG, m::Measurements, sign::Int, confs::CTInt_Confs, M::SampleMatrix; with_τshift=true)

"""


function measure_τ!(rng::AbstractRNG, m::Measurements, sign::Int, 
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
    GImp_τ = Vector{ComplexF64}(undef, length(GWeiss.data))
    for i in 1:length(GImp_τ)
        τi = GWeiss.τGrid[i]
        GImp_τ[i] = GWeiss.data[i] - 
                        τIntegrate(τi_sample -> draw_Gτ(GWeiss, τi - τi_sample),
                                    m.samples .* m.τWeights ./ m.NSamples, m.τGrid)
    end
    return GImp_τ
end

#TODO: frequency measurement
