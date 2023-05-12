@testset "solve" begin
    @test all(impSolve_IPT(GW, GW, νnGrid, 0.0, 1.0) .≈ zeros(length(νnGrid)))
end
