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
    pf = fit(xval, yval, 2) |> p -> round.(coeffs(p), digits=7) |> Polynomial
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


kG = gen_kGrid("3Dsc-0.2041241452319315", 20);

βList = Float64.(10:25:600)
UList = union(0.5:0.2:7.0)
Nit   = 15
ν0    = Nν+1 


# TESTS
U = 1.1
β = 11.1
μ = U/2 
νnGrid = jDMFT.iν_array(β, -Nν:Nν-1)
GW_data = init_weissgf_U0guess(Nν, μ, β)
GW_in = MatsubaraFunction(jDMFT.PH_transform(GW_data, U), β, νnGrid)
GW_data_2 = init_weissgf_ALguess(νnGrid, μ)
GW_in_2 = MatsubaraFunction(jDMFT.PH_transform(GW_data_2, U), β, νnGrid)


Σ_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
Σ0_Test = Array{Float64,2}(undef, length(βList), length(UList))
Σ_Test2 = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
Σ0_Test2 = Array{Float64,2}(undef, length(βList), length(UList))
Σ_AL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
Σ_NIL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
G_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
G_AL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)
G_NIL_Test = Array{ComplexF64,3}(undef, length(βList), length(UList), Nτ)

n = 1.0
ii = 0
NLoops = length(βList) * length(UList)
for (i,β) in enumerate(βList)
    νnGrid = jDMFT.iν_array(β, -Nν:Nν-1)
    for (j,U) in enumerate(UList)
        μ = U/2.0
        GW_data = init_weissgf_U0guess(Nν, μ, β)
        GW_in = MatsubaraFunction(jDMFT.PH_transform(GW_data, U), β, νnGrid)
        GW_data_2 = init_weissgf_ALguess(νnGrid, μ)
        GW_in_2 = MatsubaraFunction(jDMFT.PH_transform(GW_data_2, U), β, νnGrid)
        GW_n, Σ  = DMFTLoop(Nit, GW_in, GW_in, kG, U, n, 0.0, β * 1.0);
        GW_n_2, Σ_2  = DMFTLoop(Nit, GW_in_2, GW_in_2, kG, U, n, 0.0, β * 1.0);
        Σ_Test[i,j,:] = Σ
        Σ0_Test[i,j] = fit_Σ0(νnGrid, Σ)
        Σ_Test2[i,j,:] = Σ_2
        Σ0_Test2[i,j] = fit_Σ0(νnGrid, Σ_2)
    # println("U/β = ", U, "/", β, " // Σ0 = $(round(Σ0_Test[i,j], digits=4)) -- Σ0_2 = $(round(Σ0_Test2[i,j], digits=4))", " // 1/ν [expected ?=? fit]: ", round(U^2 * 0.5 * (1 - 0.5),digits=4), " ?=? " ,round(-((imag(Σ) .* imag(νnGrid))[end]), digits=4))
        #GW_NIL, G_NIL, Σ_NIL = solve_NonIntLimit(νnGrid, kG, μ, n, β, U)
        #GW_AL,  G_AL,  Σ_AL  = solve_AtomicLimit(νnGrid, μ, n, β, U)
        #Σ_AL_Test[i,j,:] = Σ_AL
        #Σ_NIL_Test[i,j,:] = Σ_NIL
        #G_AL_Test[i,j,:] = GW_AL
        #G_NIL_Test[i,j,:] = GW_NIL
        #G_Test[i,j,:] = GW_n.data
        global ii = ii+1
        print("\rCalculating: $(rpad(lpad(string(round(100 * ii/NLoops ,digits=2)),2),8)) % done!")
    end
end
