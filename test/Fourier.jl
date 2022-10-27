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

@testset "FT Tests" begin
    τGrid, τWeights = jDMFT.riemann(0,β-1/length(G_τ.data),length(G_τ.data))
    G_τ = ω_to_τ(GW, τWeights, τGrid)
    G2 = τ_to_ω(G_τ, νnGrid)
    @test maximum(abs.((GW.data .- G2.data) ./ GW.data)) < 1e-8
    @test real(G2.data[1] .* νnGrid[1]) ≈ 1.0
end
