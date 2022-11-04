using Pkg
Pkg.activate(".")
using jDMFT
using FastGaussQuadrature, DelimitedFiles
using Dispersions
using Polynomials 

include("test/helper_functions.jl")

function fit_Σ0(νnGrid, Σ)
    ν0_i = findmin(abs.(νnGrid))[2] + 1
    xval = imag(νnGrid[ν0_i:(ν0_i+3)])
    yval = imag(Σ[ν0_i:(ν0_i+3)])
    pf = fit(xval, yval, 2) |> p -> round.(coeffs(p), digits=4) |> Polynomial
    #println("best fit: ", pf)

    return pf.coeffs[1]
end

t_fac = 1.0 #2*sqrt(6)
# U = 1.0 #0.20412414523193154 #1.1
# μ = 0.5 #0.10206207261596577 #0.55
# β  = 10.0 #0.20412414523193154 #10.0
n  = 1.0
Nν = 200
Nτ = 400


kG = gen_kGrid("3Dsc-1.0", 20);

βList = Float64.(1:5:300)
UList = union(0.1:0.2:20.0)
Nit   = 30
ν0    = Nν+1 

ΣIm_0 = Matrix{Float64}(undef, length(βList), length(UList))
ΣIm_AL_0 = Matrix{Float64}(undef, length(βList), length(UList))
#tt = DMFTLoop(20, GW, GW, kG, U, n, μ, β*1.0)



Σ_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
Σ0_Test = Array{Float64,2}(undef, length(βList), length(UList))
Σ_AL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
Σ_NIL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
G_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
G_AL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
G_NIL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)

n = 1.0
NLoops = length(βList) * length(UList)
for (i,β) in enumerate(βList)
    νnGrid = jDMFT.iν_array(β, -Nν:Nν-1)
    for (j,U) in enumerate(UList)
        μ = U/2.0
        GW_data = init_weissgf_U0guess(Nν, μ, β)
        GW_in = MatsubaraFunction(jDMFT.PH_transform(GW_data, U), β, νnGrid)
        GW_NIL, G_NIL, Σ_NIL = solve_NonIntLimit(νnGrid, kG, μ, n, β, U)
        GW_AL,  G_AL,  Σ_AL  = solve_AtomicLimit(νnGrid, μ, n, β, U)
        GW_n, Σ  = DMFTLoop(Nit, GW_in, GW_in, kG, U, n, 0.0, β * 1.0);
        Σ_Test[i,j,:] = Σ
        Σ0_Test[i,j] = fit_Σ0(νnGrid, Σ)
    println("U/β = ", U, "/", β, " // Σ0 = $(round(Σ0_Test[i,j], digits=4))", " // 1/ν [expected ?=? fit]: ", round(U^2 * 0.5 * (1 - 0.5),digits=4), " ?=? " ,round(-((imag(Σ) .* imag(νnGrid))[end]), digits=4))
        Σ_AL_Test[i,j,:] = Σ_AL
        Σ_NIL_Test[i,j,:] = Σ_NIL
        G_AL_Test[i,j,:] = GW_AL
        G_NIL_Test[i,j,:] = GW_NIL
        G_Test[i,j,:] = GW_n.data
    end
end

# for (i,β) in enumerate(βList)
#     for (j,U) in enumerate(UList)
#         μ = U/2 
#         Σ, ΣAL = DMFTLoop(Nit, GW, GW, kG, U, n, μ, β*1.0)
#         ΣIm_0[i,j], ΣIm_AL_0[i,j] = -imag(Σ[ν0]),  -imag(ΣAL[ν0])
#     end
# end
