@testset "τIntegrate" begin
    τGrid_red, τWeights = FastGaussQuadrature.gaussradau(1000);
    τGrid = τGrid_red .* ((1.7 - 1.5) / 2) .+  ((1.7 + 1.5) / 2) 
    @test all(τGrid .== jDMFT.τGridTransform(1.5,1.7, τGrid_red))
    data = τGrid .^ 2
    fτ = τFunction(data, 1.0, τGrid, τWeights)
    @test isapprox(τIntegrate(x->x^2, τWeights, τGrid), 0.51266, rtol=1e-4)
    @test isapprox(τIntegrate(x->x^3, τWeights, τGrid), 0.8224, rtol=1e-4)
    @test isapprox(τIntegrate(x->x^2, τWeights, τGrid), τIntegrate(fτ), atol=1e-9)
    @test isapprox(τIntegrate(τGrid .^ 2, τWeights, τGrid), τIntegrate(fτ), atol=1e-9)
    rG, rW = jDMFT.riemann(3.13, 5.67, 1000)
    grG, grW = jDMFT.gaussradau(3.13, 5.67, 1000)

    @test τIntegrate(x-> x^5, rW, rG) ≈ 5381.21 rtol=0.1
    @test τIntegrate(x-> x^5, grW, grG) ≈ 5381.21 rtol=0.1
end

@testset "naive TF" begin
    τGrid_riemann,τWeights_riemann = jDMFT.riemann(0,β-0.0001,length(G_τ.data))
    G_τ_riemann  = ω_to_τ(GW, ones(Float64, length(τGrid_riemann)) ./ length(τGrid_riemann), τGrid_riemann)
    G2_riemann = τ_to_ω(G_τ_riemann, νnGrid, [1.0], [-0.5])
    G2_GR = τ_to_ω(G_τ, νnGrid, [1.0], [-0.5])
    @test maximum(abs.((GW.data .- G2_riemann.data) ./ GW.data)) < 1.0
    @test maximum(abs.((GW.data .- G2_GR.data) ./ GW.data)) < 1/(length(G_τ.data)^1.8)
end
