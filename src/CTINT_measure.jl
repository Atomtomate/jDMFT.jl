#TODO: use more advanced data structure to capture statistics
mutable struct Measurements
    samples::Vector{ComplexF64}
    τWeights::Vector{Float64}
    τGrid::Vector{Float64}
    NSamples::Int
    totalSign::Int
end

"""
    CTInt_Confs


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


#TODO: test and implement fast-version (only draw random indices)
#TODO: replace linear interpolation by some cubic lin. inter, using legndre grid 
function draw_Gτ(GWeiss::τFunction, τi::Float64)::ComplexF64
    Nτ = length(GWeiss.data)
    sign = τi < 0 ? -1 : 1
    τi = τi < 0 ? τi + GWeiss.β  : τi
    ind = searchsortedfirst(GWeiss.τGrid, τi)
    ind = ind > Nτ ? ind-1 : ind
    m = if ind < Nτ 
        (GWeiss.data[ind+1] - GWeiss.data[ind])/(GWeiss.τGrid[ind+1] - GWeiss.τGrid[ind]) 
    else
        (GWeiss.data[1] - GWeiss.data[ind])/(GWeiss.τGrid[1] - GWeiss.τGrid[ind] + GWeiss.β)
    end
    val = GWeiss.data[ind] + m * (τi - GWeiss.τGrid[ind])
    return val
end


function measure_τ!(rng::AbstractRNG, m::Measurements, success::Bool, sign::Int, 
                  C::CTInt_Confs, M::SampleMatrix)
    for i in 1:M.N
        τi = C.τiList[i] - τ_shift
        M.rowCache[i] = draw_Gτ(C.τList, C.GWeiss, τi, C.β)
    end
    M.colCache[:] = view(M.data, 1:M.N, 1:M.N) * view(M.rowCache, 1:M.N)
    N_bins = length(m.samples)
    for i in 1:M.N
        τ_shift = rand(rng, Float64) * C.β
        τi = C.τList[i] - τ_shift
        s = 2*(τi > 0) - 1
        τi  = τi + (τi < 0)*C.β
        ind = ceil(Int, τi/C.β * N_bins)
        m.samples[ind] += sign * s * M.colCache[i]
    end
    m.NSamples += 1
    m.totalSign += sign
end

# struct τFunction    
#    14     data::Vector{ComplexF64}
#    13     β::Float64
#    12     τGrid::AbstractVector
#    11     τWeights::AbstractVector
#    10 end
function measure_GImp_τ(m::Measurements, GWeiss::τFunction, Niw)
    data = Vector{ComplexF64}(undef, length(GWeiss.data))
    for i in 1:length(data)
        g0_data = [draw_Gτ(GWeiss) for τi in m.τGrid]
        for j in 1:length(m.samples)

        end
    end
    #TODO: fourier trafo for each iν
end


#TODO: frequency measurement
#
#
