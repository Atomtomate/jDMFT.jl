@testset "τIndex" begin
    @test all(jDMFT.τIndex(2.0, [1.0,2.1,3.0], 3.1) .== (1,2))
    @test all(jDMFT.τIndex(3.2, [1.0,2.1,3.0], 3.1) .== (-1,1))
    @test all(jDMFT.τIndex(0.02, [1.0,2.1,3.0], 3.1) .== (1,1))
    @test all(jDMFT.τIndex(0.00, [1.0,2.1,3.0], 3.1) .== (1,1))
    @test all(jDMFT.τIndex(-0.01, [1.0,2.1,3.0], 3.1) .== (-1,3))
end

@testset "draw_Gτ" begin
    ii = collect(0:0.01:(2*π-0.01));
    v = sin.(ii);
    G0_test = τFunction(convert.(ComplexF64,v), 2π, ii, 1/length(ii) .* ones(length(ii)))
    @test jDMFT.draw_Gτ(G0_test, 0.1) == v[11]
    @test real(jDMFT.draw_Gτ(G0_test, 0.105)) > v[11] && real(jDMFT.draw_Gτ(G0_test, 0.105)) < v[12]
    @test jDMFT.draw_Gτ(G0_test, 2π) ≈ 0
    @test jDMFT.draw_Gτ(G0_test, -1.0) ≈ - jDMFT.draw_Gτ(G0_test, -1.0+2π)
    @test jDMFT.draw_Gτ(G0_test, -1.0 - 2*2π) ≈ jDMFT.draw_Gτ(G0_test, -1.0+2*2π)
end
