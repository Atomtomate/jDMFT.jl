@testset "CTINT measure" begin
    include("CTINT_measure.jl")
end

@testset "CTInt_Confs" begin
    G0_ind = [0, 0.01, 0.02, 2.2, 5.5, 10.0, 15.5, 20.0] 
    G0 = τFunction(convert.(ComplexF64,collect(1:length(G0_ind))), 20.0, G0_ind, ones(length(G0_ind))/length(G0_ind))
    c = jDMFT.CTInt_Confs(G0, G0_ind, 1.1)
    @test c.β ≈ 20.0
    @test c.U ≈ 1.1
    @test length(c.τList) == length(G0_ind)
    @test length(c.τiList) == 10
    @test length(c.siList) == 10
end

@testset "update cache" begin
    G0_ind = [0, 0.01, 0.02, 2.2, 5.5, 10.0, 15.5, 20.0] 
    G0 = τFunction(convert.(ComplexF64,collect(1:length(G0_ind))), 20.0, G0_ind, ones(length(G0_ind))/length(G0_ind))
    c = jDMFT.CTInt_Confs(G0, G0_ind, 1.1)
    M = SampleMatrix(15)
    @test jDMFT.update_cache!(c,M) == 0
    @test length(M.rowCache) == 15
    M.N = 15
    @test jDMFT.update_cache!(c,M) == 1
    @test length(M.rowCache) == 25
    @test length(M.colCache) == 25
    @test all(size(M.data) .== (25,25))
    @test length(c.τiList) == 25
    @test length(c.siList) == 25
end

@testset "try_insert!" begin
    rng = MersenneTwister(0)
    cf = jDMFT.CTInt_Confs(G_τ, U)
    M = SampleMatrix(15)
end

@testset "sample" begin
end


@testset "sample_step" begin
end
