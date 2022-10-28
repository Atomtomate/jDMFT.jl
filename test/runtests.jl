using Test, Logging
using Random
using FastGaussQuadrature
# using Dispersions
using jDMFT
using DelimitedFiles

debuglogger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debuglogger)

## Example ED hybridization values for U=1, β=20.6
include("helper_functions.jl")

# r = readdlm("test_data/g0iw.txt")
# GWeiss_w2dyn = r[3:2002,6] .+ r[2003:end,6]) ./2 .+ 1im .* (r[3:2002,7] .+ r[2003:end,7])
# iw_w2dyn = r[802:1200,5] .* 1im

t_fac = 1#2*sqrt(6)
U = 2.0 * t_fac
β  = 14.0 * t_fac
μ = 1.00 * t_fac
n = 1.0
# U = 0.7 * t_fac
# β  = 1.0 * t_fac
# μ = 0.35 * t_fac
Nν = 200
Nτ = 400
νnGrid = jDMFT.iν_array(β, -Nν:Nν-1)
#τGrid_red, τWeights = gaussradau(Nτ);
#τGrid = τGrid_red .* β ./ 2 .+ β ./ 2;
GW_in = init_weissgf_from_ed((@__DIR__)*"/test_data/U2.0_b14.0_hubb.andpar", Nν)
# GW_in = init_weissgf_from_ed((@__DIR__)*"/test_data/U0.7_b1.0_hubb.andpar", Nν)
GW = MatsubaraFunction(jDMFT.PH_transform(GW_in, U), β, νnGrid, [0.0, 1.0])
τGrid, τWeights = jDMFT.riemann(0.0,β-1/Nτ,Nτ)
G_τ = ω_to_τ(GW, τWeights, τGrid)
#GImp_test = init_GImp((@__DIR__)*"/test_data/U2.0_b14.0_gm_wim", Nν)
#GImp_τ_test = ω_to_τ(GW, τWeights, τGrid)

include("Random.jl")

@testset "DataTypes" begin
    include("DataTypes.jl")
end

@testset "Fourier" begin
    include("Fourier.jl")
    include("FourierFFT.jl")
end

@testset "GFTools" begin
    include("GFTools.jl")
end

@testset "DMFT Loop" begin
    include("DMFTLoop.jl")
end

@testset "IPT" begin
    include("IPT.jl")
end

@testset "CTINT" begin
    include("CTINT.jl")
end
