include("CTINT_measure.jl")

#TODO: use some sort of module instead...
function update_cache!(C::CTInt_Confs, M::SampleMatrix; threshold = 0, incr = 10)
    if M.N - length(M.rowCache) >= threshold
        lNew = size(M.data,1) + incr
        data_cache = similar(M.data, (size(M.data) .+ incr)...)
        data_cache[axes(M.data)...] .= M.data
        M.data = data_cache
        resize!(M.rowCache,lNew)
        resize!(M.colCache,lNew)
        resize!(C.τiList,lNew)
        resize!(C.siList,lNew)
        return 1
    else
        return 0
    end
end

function sample!(measurements::Measurements, N::Int, β::Float64, U::Float64, GWeiss::Vector, τGrid::Vector)
    rng = MersenneTwister(0)
    confs = CTInt_Confs(β, U, GWeiss, τGrid, 15)
    M = SampleMatrix(15)
    for i in 1:N
        success, sign = sample_step!(rng, confs, M, GWeiss)
        measure_τ!(rng, measurements, success, sign, C, M)
    end
end

function sample_step!(rng::T, C::CTInt_Confs, M::SampleMatrix) where  T  <: AbstractRNG
    step = rand(rng, [:insert, :remove])
    sign = 0
    success = false
    if step == :insert
        success, sign = try_insert!(rng, C, M)
    elseif step == :remove
        success, sign = try_remove!(rng, C, M)
    end
    return success, sign
end

function try_insert!(rng::AbstractRNG, C::CTInt_Confs, M::SampleMatrix; δ::Float64 = 0.0001)
    update_cache!(C, M)
    #TODO: movethis to config, do not force paramagnetic solution / G_σ = G_-σ
    σ = 1
    N = M.N
    Nnew = M.N + 1
    τi = rand(rng, Float64) * C.β
    si = rand(rng, [0,1])
    ζi = rand(rng, Float64)
#    @debug τi si ζi
    M.S = C.GWeiss[1] - 0.5 + σ*si*(0.5 + δ)
    for j in 1:M.N
        τj = C.τiList[j]
        M.rowCache[j] = draw_Gτ(C.GWeiss, τi - τj)
        M.colCache[j] = draw_Gτ(C.GWeiss, τj - τi)
    end
    Mi = view(M.data, 1:N, 1:N)
    Ri = transpose(view(M.rowCache, 1:N))
    Ci = view(M.colCache, 1:N)
    S = N > 0 ? 1/(M.S .- Ri*(Mi * Ci)) : M.S
    R = real((1/S) * (1/S) *(-C.β*C.U/(Nnew)))
#    @debug Mi Ci S R
    success = false
    if ζi < abs(R) # Update successfull
        M.data[1:N,Nnew]  = -(Mi * Ci) * M.S
        M.data[Nnew,1:N]  = -M.S * (Ri * Mi)
        M.data[1:N,1:N]   = Mi + (Mi*Ci)*M.S*(Ri*Mi)
        M.data[Nnew,Nnew] = S
        #TODO update M
        M.N = Nnew
        C.τiList[Nnew] = τi
        C.siList[Nnew] = si
        success = true
    else
        #TODO: Accumulate with shifted tau
    end
    return success, R
end


function try_remove!(rng::AbstractRNG, C::CTInt_Confs, M::SampleMatrix)
    success = false
    N = M.N
    S = M.data[N,N]
    Ri = transpose(view(M.rowCache, 1:(N-1)))
    Ci = view(M.colCache, 1:(N-1))
    ζi = rand(rng, Float64)
    R = (S*S)/(-C.β*C.U/(M.N))
    if ζi < abs(R) 
        N = N - 1
        M.data[1:N,1:N] = M.data[1:N,1:N] .- (Ci .* Ri) ./ S
        M.N = N
    else
        #TODO: Accumulate with shifted tau
    end
    return success, R
end
