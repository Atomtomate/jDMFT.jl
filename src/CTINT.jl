# ==================================================================================================== #
#                                             CTINT.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Julian Stobbe                                                                    #
#   Last Edit Date  : 10.06.22                                                                         #
# ----------------------------------------- Description ---------------------------------------------- #
#   Interaction expansion impurity solver.                                                             #
# -------------------------------------------- TODO -------------------------------------------------- #
#   -  Up/Down SampleMatrix are hardcoded for now!                                                     #
# ==================================================================================================== #


include("CTINT_measure.jl")

#TODO: hardcoded GWeiss_up == GWeiss_do
function setup_CTINT(GWeiss::τFunction, sample_τGrid::Vector{Float64}, U::Float64)
    cf = CTInt_Confs(GWeiss, U)
    #TODO: hardcoded up/down
    me = [Measurements(zeros(ComplexF64, length(sample_τGrid)), sample_τGrid, 0, 0, 0), 
          Measurements(zeros(ComplexF64, length(sample_τGrid)), sample_τGrid, 0, 0, 0)] 
    MList = [SampleMatrix(), SampleMatrix()]
    return cf, me, MList
end

function sample(GWeiss::τFunction, U::Float64, N::Int; N_warmup::Int=0, sample_τGrid=nothing, with_τshift=true)
    sample_τGrid = sample_τGrid == nothing ? GWeiss.τGrid : sample_τGrid
    cf, me, MList = setup_CTINT(GWeiss, sample_τGrid, U)
    warmup!(cf, MList, N_warmup)
    GImp_test, GImp_test2 = sample!(cf, me, MList, N, with_τshift=with_τshift)
    return cf, me, MList, GImp_test, GImp_test2
end

function warmup!(confs::CTInt_Confs, MList::Vector{SampleMatrix}, N::Int)
    rng = MersenneTwister(0)
    for i in 1:N
        success, sign = sample_step!(rng, confs, MList)
    end
end

function sample!(confs::CTInt_Confs, measurements::Vector{Measurements}, MList::Vector{SampleMatrix}, N::Int; with_τshift=true)
    rng = MersenneTwister(1)

    GImp_τTest = [zeros(ComplexF64, length(confs.GWeiss.data)), zeros(ComplexF64, length(confs.GWeiss.data))]
    GImp_νTest = [zeros(ComplexF64, 200), zeros(ComplexF64, 200)]
    νnGrid = iν_array(confs.β, -100:99)
    totalSign = 0
    for i in 1:N
        success, sign = sample_step!(rng, confs, MList)
        totalSign += sign
        # if sign == -1
        #     println("WARNING: something went wrong on i=$i")
        #     println("M[1] = ", MList[1].N)
        #     println("M[2] = ", MList[2].N)
        #     println("cconfs = ", confs.τiList[1:MList[1].N])
        #     println("cconfs = ", confs.siList[1:MList[1].N])
        # end
        #println("loop:", sign)
        for i in 1:length(MList)
            accumulate!(rng, measurements[i], sign, confs, MList[i], with_τshift=with_τshift)
            measure_GImp_ωn!(GImp_νTest[i], νnGrid, MList[i], confs, sign)
            #measure_GImp_τ!(GImp_τTest[i], MList[i], confs.GWeiss, sign)
        end
    end
    for i in 1:length(MList)
        GImp_τTest = GImp_τTest ./ totalSign
        GImp_νTest = GImp_νTest ./ totalSign 
    end
    return GImp_τTest, GImp_νTest
end

function sample_step!(rng::T, confs::CTInt_Confs, MList::Vector{SampleMatrix}) where  T  <: AbstractRNG
    step = rand(rng, [:insert, :remove])
    sign = 0
    success = false
    if step == :insert
        success, sign = try_insert!(rng, confs, MList)
    elseif step == :remove
        success, sign = try_remove!(rng, confs, MList)
    end
    return success, sign
end

function try_insert!(rng::AbstractRNG, confs::CTInt_Confs, MList::Vector{SampleMatrix})
    success = false
    τi = rand(rng, Float64) * confs.β
    si = rand(rng, [-1,1])
    ζi = rand(rng, Float64)
    #TODO: move this to config, do not force paramagnetic solution / G_σ = G_-σ
    σList = [-1,1]

    N = MList[1].N
    R = -confs.β * confs.U/(N+1)
    for (σ,M) in zip(σList, MList)
        update_cache!(confs, M)
        # G
        M.S = confs.GWeiss.data[1] -  (0.5 + σ*si*(0.5+confs.δ)) #(0.25 + σ*si*(0.25+confs.δ))
        for j in 1:N
            τj = confs.τiList[j]
            M.rowCache[j] = draw_Gτ(confs.GWeiss, τi - τj) #TODO: draw from up/down GF, after generalization
            M.colCache[j] = draw_Gτ(confs.GWeiss, τj - τi)
        end
        Mi = view(M.data, 1:N, 1:N)
        Ri = transpose(view(M.rowCache, 1:N))
        Ci = view(M.colCache, 1:N)
        S = N > 0 ? 1/(M.S - Ri*(Mi * Ci)) : 1/M.S
        # TODO: this is computed again in fast_update
        R = R * (1/S)
        if DEBUG
            confs.τiList[N+1] = τi
            confs.siList[N+1] = si
            pd = plain_det(N, N+1, confs)
            println("insert (N=$N): plain det: $(round.(real(pd),digits=6)) vs. fast $(round(real(1/S),digits=6))")
        end
    end

    #println("ins: ", 2*(real(R)), " // ", R)
    if ζi < abs(real(R)) # Update successfull
        for M in MList
            fast_update_incr!(M, N+1)
        end
        confs.τiList[N+1] = τi
        confs.siList[N+1] = si
        success = true
    end
    return success, 2*(real(R) > 0) - 1
end


function try_remove!(rng::AbstractRNG, confs::CTInt_Confs, MList::Vector{SampleMatrix})
    N = MList[1].N
    success = false
    if N == 0
        return false, 1 # TODO: fixed sign to -1 for G(0^+)// +1 for G(0^-) may lead to bugs later
    end

    R = prod(M.data[N,N] for M in MList)
    ζi = rand(rng, Float64)

    R = -N*R/(confs.β*confs.U)
    if DEBUG
        pd = 1 / plain_det(N-1, N, confs)
        println("remove (N=$N): plain det: $(round.(real(pd),digits=6)) vs. fast $(round(real(S),digits=6))")
    end
    if ζi < abs(real(R))
        for M in MList
            fast_update_decr!(M, N-1)
        end
        success = true
    end
    return success, 2*(real(R) > 0) - 1
end



# ================================================================================
#                           Testing and stability                                =
# ================================================================================

function rebuild_SampleMatrix(N::Int, confs::CTInt_Confs, σ)
    new_data = Matrix{ComplexF64}(undef, N, N)
    for i in 1:N
        τi = confs.τiList[i]
        si = confs.siList[i]
        for j in 1:N
            τj = confs.τiList[j]
            new_data[i,j] = draw_Gτ(confs.GWeiss, τi - τj)
        end
        new_data[i,i] -= (0.5 + σ*si*(0.5+confs.δ))
    end
    return inv(new_data)
end

function plain_det(Nold::Int, Nnew::Int, confs::CTInt_Confs)
    Mold = rebuild_SampleMatrix(Nold, confs)
    Mnew = rebuild_SampleMatrix(Nnew, confs)
    return det(inv(Mnew))/det(inv(Mold))
end
