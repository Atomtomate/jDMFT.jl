@testset "τIntegrate" begin
    τGrid_red, τWeights = gaussradau(1000);
    τGrid = τGrid_red .* ((1.7 + 1.5) / 2) .+  ((1.7 - 1.5) / 2) 
    data = τGrid .^ 2
    fτ = τFunction(data, 1.0, τGrid, τWeights)
    @test isapprox(τIntegrate(x->x^2, τWeights, τGrid), 2.76267, rtol=1e-4)
    @test isapprox(τIntegrate(x->x^3, τWeights, τGrid), 0.8224, rtol=1e-4)
    @test isapprox(τIntegrate(x->x^2, τWeights, τGrid), τIntegrate(fτ), atol=1e-9)
end

@testset "naive TF" begin
    τGrid_riemann = collect(LinRange(0,β-0.0001, length(G_τ.data)))
    G_τ_riemann  = ω_to_τ(GW, ones(Float64, length(τGrid_riemann)) ./ length(τGrid_riemann), τGrid_riemann)
    G2_riemann = τ_to_ω(G_τ_riemann, νnGrid, [1.0], [-0.5])
    G2_GR = τ_to_ω(G_τ, νnGrid, [1.0], [-0.5])
    @test maximum(abs.((GW.data .- G2_riemann.data) ./ GW.data)) < 1.0
    @test maximum(abs.((GW.data .- G2_GR.data) ./ GW.data)) < 1e-6
end
