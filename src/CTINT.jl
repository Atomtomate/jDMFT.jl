include("CTINT_measure.jl")

function setup_CTINT(GWeiss::τFunction, sample_τGrid::Vector{Float64}, U::Float64)
    cf = CTInt_Confs(GWeiss, U)
    me = Measurements(zeros(ComplexF64, length(sample_τGrid)), sample_τGrid, 0, 0)
    M = SampleMatrix()
    return cf, me, M
end

function sample(GWeiss::τFunction, U::Float64, N::Int; sample_τGrid=nothing)
    sample_τGrid = sample_τGrid == nothing ? GWeiss.τGrid : sample_τGrid
    cf, me, M = setup_CTINT(GWeiss, sample_τGrid, U)
    sample!(cf, me, M, N)
end

function sample!(confs::CTInt_Confs, measurements::Measurements, M::SampleMatrix, N::Int)
    rng = MersenneTwister(0)
    for i in 1:N
        success, sign = sample_step!(rng, confs, M)
        accumulate!(rng, measurements, sign, confs, M)
    end
end

function sample_step!(rng::T, confs::CTInt_Confs, M::SampleMatrix) where  T  <: AbstractRNG
    step = rand(rng, [:insert, :remove])
    sign = 0
    success = false
    if step == :insert
        success, sign = try_insert!(rng, confs, M)
    elseif step == :remove
        success, sign = try_remove!(rng, confs, M)
    end
    return success, sign
end

function try_insert!(rng::AbstractRNG, confs::CTInt_Confs, M::SampleMatrix; δ::Float64 = 0.0001)
    update_cache!(confs, M)
    #TODO: movethis to config, do not force paramagnetic solution / G_σ = G_-σ
    σ = 1
    N = M.N
    Nnew = M.N + 1
    τi = rand(rng, Float64) * confs.β
    si = rand(rng, [0,1])
    ζi = rand(rng, Float64)
#    @debug τi si ζi
    M.S = confs.GWeiss.data[1] - 0.5 + σ*si*(0.5 + δ)
    for j in 1:M.N
        τj = confs.τiList[j]
        M.rowCache[j] = draw_Gτ(confs.GWeiss, τi - τj)
        M.colCache[j] = draw_Gτ(confs.GWeiss, τj - τi)
    end
    Mi = view(M.data, 1:N, 1:N)
    Ri = transpose(view(M.rowCache, 1:N))
    Ci = view(M.colCache, 1:N)
    S = N > 0 ? 1/(M.S .- Ri*(Mi * Ci)) : M.S
    R = real((1/S) * (1/S) *(-confs.β*confs.U/(Nnew)))
#    @debug Mi Ci S R
    success = false
    if ζi < abs(R) # Update successfull
        M.data[1:N,Nnew]  = -(Mi * Ci) * M.S
        M.data[Nnew,1:N]  = -M.S * (Ri * Mi)
        M.data[1:N,1:N]   = Mi + (Mi*Ci)*M.S*(Ri*Mi)
        M.data[Nnew,Nnew] = S
        #TODO update M
        M.N = Nnew
        confs.τiList[Nnew] = τi
        confs.siList[Nnew] = si
        success = true
    else
        #TODO: Accumulate with shifted tau
    end
    return success, 2*(real(R) > 0) - 1
end


function try_remove!(rng::AbstractRNG, confs::CTInt_Confs, M::SampleMatrix)
    success = false
    N = M.N
    S = M.data[N,N]
    Ri = transpose(view(M.rowCache, 1:(N-1)))
    Ci = view(M.colCache, 1:(N-1))
    ζi = rand(rng, Float64)
    R = (S*S)/(-confs.β*confs.U/(M.N))
    if ζi < abs(R) 
        N = N - 1
        M.data[1:N,1:N] = M.data[1:N,1:N] .- (Ci .* Ri) ./ S
        M.N = N
    else
        #TODO: Accumulate with shifted tau
    end
    return success, 2*(real(R) > 0) - 1
end
