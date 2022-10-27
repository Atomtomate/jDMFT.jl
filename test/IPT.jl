@testset "solve" begin
    t1 = τFunction(G_τ.data .^ 3, β, τGrid, τWeights, [0.5^3, 0.0])
    t1_ft = τ_to_ω(t1, νnGrid);
    @test all(impSolve_IPT(G_τ, G_τ, νnGrid, 0.0, 1.0) .≈ zeros(length(νnGrid)))
end
