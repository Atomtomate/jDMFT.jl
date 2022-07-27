@testset "update cache" begin
    G0_ind = [0, 0.01, 0.02, 2.2, 5.5, 10.0, 15.5, 20.0] 
    G0 = τFunction(convert.(ComplexF64,collect(1:length(G0_ind))), 20.0, G0_ind, ones(length(G0_ind))/length(G0_ind))
    c = jDMFT.CTInt_Confs(G0, G0_ind, 1.1)
    M = jDMFT.SampleMatrix(15)
    @test jDMFT.update_cache!(c,M) == 0
    @test length(M.rowCache) == 15
    M.N = 15
    M.data[:,:] = randn(ComplexF64, size(M.data)...)
    old_data = deepcopy(M.data)
    old_row_data = deepcopy(M.rowCache)
    old_col_data = deepcopy(M.colCache)
    @test jDMFT.update_cache!(c,M) == 1
    @test length(M.rowCache) == 25
    @test length(M.colCache) == 25
    @test all(size(M.data) .== (25,25))
    @test length(c.τiList) == 25
    @test length(c.siList) == 25
    @test all(M.data[1:15,1:15] .≈ old_data)
    @test all(M.rowCache[1:15] .≈ old_row_data)
    @test all(M.colCache[1:15] .≈ old_col_data)
end


@testset "Measurements" begin
    β_test = 2.2
    τGrid_red, τWeights = gaussradau(10);
    τGrid = τGrid_red .* β_test ./ 2 .+ β_test ./ 2;
    τGrid_riemann = collect(LinRange(0,β_test-0.001/10, 10))

    sample_τGrid = 1.1 .* collect(1:10)
    sample_τWeights = 1.1 .* collect(6:15)
    m = jDMFT.Measurements(zeros(ComplexF64, length(sample_τGrid)), sample_τWeights, sample_τGrid, 0, 0)
    @test all(m.samples .== zeros(ComplexF64, length(sample_τGrid)))
    @test all(m.τWeights .≈ sample_τWeights)
    @test all(m.τGrid .≈ sample_τGrid)
    @test m.NSamples == 0
    @test m.totalSign == 0
    m2 = jDMFT.Measurements(10, β_test, :GaussQuad)
    @test m2.NSamples == 0
    @test m2.totalSign == 0
    @test sum(m2.samples) == 0
    @test all(m2.τGrid .≈ τGrid)
    @test all(m2.τWeights .≈ τWeights)
    m3 = jDMFT.Measurements(10, β_test, :Riemann)
    @test m3.NSamples == 0
    @test m3.totalSign == 0
    @test sum(m3.samples) == 0
    @test all(m3.τGrid .≈ τGrid_riemann)
    @test all(m3.τWeights .≈ ones(Float64, length(10)) ./ length(10))
    @test_throws DomainError jDMFT.Measurements(10, β_test, :throw) 
end

@testset "measure_τ" begin
    rng = MersenneTwister(0)
    cf = jDMFT.CTInt_Confs(G_τ, U)
    N = length(cf.τiList)
    m_R = jDMFT.Measurements(N, G_τ.β, :Riemann)
    m_GQ = jDMFT.Measurements(N, G_τ.β, :GaussQuad)
    M = jDMFT.SampleMatrix()
    jDMFT.measure_τ!(rng, m_R, -1, cf, M, with_τshift=false)
    jDMFT.measure_τ!(rng, m_GQ, -1, cf, M, with_τshift=false)
    @test m_R.NSamples == 1
    @test m_R.totalSign == -1
    @test sum(m_R.samples) == 0
    @test m_R.NSamples == m_GQ.NSamples
    @test m_R.totalSign == m_GQ.totalSign
    @test sum(m_R.samples) == sum(m_GQ.samples)
    M.N = 1
    M.data[1,1] = 1.0
    τi_test1 = 0.0
    cf.τiList[1] = τi_test1
    jDMFT.measure_τ!(rng, m_R, -1, cf, M, with_τshift=false)
    jDMFT.measure_τ!(rng, m_GQ, -1, cf, M, with_τshift=false)
    @test sum(m_R.samples) ≈ -jDMFT.draw_Gτ(G_τ, τi_test1)
    @test sum(m_R.samples) == sum(m_GQ.samples)

    m_R = jDMFT.Measurements(N, G_τ.β, :Riemann)
    m_GQ = jDMFT.Measurements(N, G_τ.β, :GaussQuad)
    M.data[1,1] = 1/G_τ.data[1]
    jDMFT.measure_τ!(rng, m_R, 1, cf, M, with_τshift=false)
    jDMFT.measure_τ!(rng, m_GQ, 1, cf, M, with_τshift=false)
    @test m_R.samples[1] ≈ 1.0
    @test m_GQ.samples[1] ≈ 1.0
    @test sum(m_R.samples) ≈ 1.0
    @test sum(m_GQ.samples) ≈ 1.0
end

@testset "measure_GImp_τ" begin
    rng = MersenneTwister(0)
    cf = jDMFT.CTInt_Confs(G_τ, U)
    m_R = jDMFT.Measurements(5, G_τ.β, :Riemann)
    m_GQ = jDMFT.Measurements(5, G_τ.β, :GaussQuad)
    M = jDMFT.SampleMatrix()
    M.N = 1
    τi_test1 = 0.0

    # Set measurement to 0, then GImp == GWeiss
    M.data[1,1] = 0.0
    cf.τiList[1] = τi_test1
    jDMFT.measure_τ!(rng, m_R, -1, cf, M, with_τshift=false)
    jDMFT.measure_τ!(rng, m_GQ, -1, cf, M, with_τshift=false)
    @test sum(jDMFT.measure_GImp_τ(m_R, G_τ)) ≈ sum(jDMFT.measure_GImp_τ(m_GQ, G_τ))
    @test sum(jDMFT.measure_GImp_τ(m_R, G_τ)) ≈ sum(G_τ.data)

    # Set measurement to -δ_τ0, then GImp == GWeiss + GWeiss(τ = 0)
    m_R = jDMFT.Measurements(5, G_τ.β, :Riemann)
    m_GQ = jDMFT.Measurements(5, G_τ.β, :GaussQuad)
    M.data[1,1] = 1/G_τ.data[1]
    jDMFT.measure_τ!(rng, m_R, 1, cf, M, with_τshift=false)
    jDMFT.measure_τ!(rng, m_GQ, 1, cf, M, with_τshift=false)
    t_R = jDMFT.measure_GImp_τ(m_R, G_τ)
    t_GQ = jDMFT.measure_GImp_τ(m_GQ, G_τ)
    println("DBG: ", round.(real.(t_R), digits=3))
    println("DBG: ", round.(real.(G_τ.data), digits=3))
    #println("DBG: ", t_GQ)
    @test sum(jDMFT.measure_GImp_τ(m_R, G_τ)) ≈ sum(jDMFT.measure_GImp_τ(m_GQ, G_τ))
    @test sum(jDMFT.measure_GImp_τ(m_R, G_τ)) ≈ 2*sum(G_τ.data)
    println(m_R.samples .* m_R.τWeights ./ m_R.NSamples)
    println(m_GQ.samples .* m_GQ.τWeights ./ m_GQ.NSamples)
    println(sum(G_τ.data .+ G_τ.data[1]))
    println(sum(G_τ.data .- G_τ.data[1]))
end
