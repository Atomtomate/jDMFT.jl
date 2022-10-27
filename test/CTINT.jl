@testset "CTINT measure" begin
    include("CTINT_measure.jl")
end

@testset "CTInt_Confs" begin
    G0_ind = [0, 0.01, 0.02, 2.2, 5.5, 10.0, 15.5, 20.0] 
    G0 = τFunction(convert.(ComplexF64,collect(1:length(G0_ind))), 20.0, G0_ind, ones(length(G0_ind))/length(G0_ind), [0.0, 0.0])
    c = jDMFT.CTInt_Confs(G0, G0_ind, 1.1)
    @test c.β ≈ 20.0
    @test c.U ≈ 1.1
    @test length(c.τList) == length(G0_ind)
    @test length(c.τiList) == 10
    @test length(c.siList) == 10
end

@testset "sample" begin
    cf,me,M = sample(G_τ, U, 15)
    M1_new = jDMFT.rebuild_SampleMatrix(M[1].N,cf, -1) #TODO: for now [-1,1] is hardcoded
    M2_new = jDMFT.rebuild_SampleMatrix(M[2].N,cf,  1)
    #println("took $(me.NSamples) with diff: ",sum(abs.(M.data[1:M.N,1:M.N] .- M_new)))
    @test all(M[1].data[1:M[1].N,1:M[1].N] .≈ M1_new)
    @test all(M[2].data[1:M[2].N,1:M[2].N] .≈ M2_new)
end
