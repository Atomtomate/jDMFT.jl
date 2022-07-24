@testset "draw_Gτ" begin
    ii = collect(0:0.01:(2*π-0.01));
    v = sin.(ii);
    G0_test = τFunction(convert.(ComplexF64,v), 2π, ii, 1/length(ii) .* ones(length(ii)))
    @test jDMFT.draw_Gτ(G0_test, 0.1) == v[11]
    @test real(jDMFT.draw_Gτ(G0_test, 0.105)) > v[11] && real(jDMFT.draw_Gτ(G0_test, 0.105)) < v[12]
    @test jDMFT.draw_Gτ(G0_test, 2π) ≈ 0
end
