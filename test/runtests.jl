using Test, Logging
using Random
using OffsetArrays
using FastGaussQuadrature
using Dispersions
using jDMFT
using DelimitedFiles

debuglogger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debuglogger)

## Example ED hybridization values for U=1, β=20.6
include("helper_functions.jl")
include("Random.jl")
include("DataTypes.jl")

U = 0.7
μ = 0.35
β = 1.0
Nν = 20
Nτ = 100
νnGrid = jDMFT.iν_array(β, -Nν:Nν-1)
τGrid_red, τWeights = gaussradau(Nτ);
τGrid = τGrid_red .* β ./ 2 .+ β ./ 2;
GW_in = init_weissgf_from_ed((@__DIR__)*"/test_data/U$(U)_b$(β)_hubb.andpar", Nν)
GW = MatsubaraFunction(jDMFT.PH_transform(GW_in, U), β, νnGrid, [1.0], [-0.5])
G_τ = ω_to_τ(GW, τWeights, τGrid)
GImp_test = init_GImp((@__DIR__)*"/test_data/U0.7_b1.0_gm_wim", Nν)
GImp_τ_test = ω_to_τ(GW, τWeights, τGrid)

@testset "Fourier" begin
    include("Fourier.jl")
end
@testset "GFTools" begin
    include("GFTools.jl")
end
@testset "DMFT Loop" begin
    include("DMFTLoop.jl")
end
@testset "CTINT" begin
    include("CTINT.jl")
end

