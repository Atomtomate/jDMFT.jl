@testset "MatsubaraFunction" begin
    t = MatsubaraFunction([1.0+0.0im], 1.0, [1.1+0.0im]) 
    @test all(t.data .== [1.0+0.0im])
    @test all(t.fGrid .== [1.1+0.0im])
    @test all(t.tail_coeffs .== [0.0, 1.0])
    @test t.β == 1.0
end

@testset "τFunction" begin
    t = τFunction([1.0+0.0im], 1.0, [1.1], [1.2]) 
    @test all(t.data .== [1.0+0.0im])
    @test all(t.τGrid .== [1.1])
    @test all(t.τWeights .== [1.2])
    @test all(t.tail_coeffs .== [0.5, 0.0])
    @test t.β == 1.0
end

@testset "tailTransform_ω_to_τ" begin
    @test all(jDMFT.tailTransform_ω_to_τ([0.0, 1.0]) .== [0.5, 0.0])
    @test_throws "Only length 2" jDMFT.tailTransform_ω_to_τ([0.0]) 
    @test_throws "constant coefficient " jDMFT.tailTransform_ω_to_τ([1.0, 0.0]) 
end

@testset "tailTransform_τ_to_ω" begin
    @test_throws "Only length 2" jDMFT.tailTransform_τ_to_ω([0.0]) 
    @test_throws "1/τ coefficient " jDMFT.tailTransform_τ_to_ω([0.0, 1.0]) 
    @test all(jDMFT.tailTransform_τ_to_ω([0.5, 0.0]) .== [0.0, 1.0])
    @test all(jDMFT.tailTransform_τ_to_ω(jDMFT.tailTransform_ω_to_τ([0.0, 1.0])) .== [0.0, 1.0])
end

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
    G0_test = τFunction(convert.(ComplexF64,v), 2π, ii, 1/length(ii) .* ones(length(ii)), [0.5, 0.0])
    @test jDMFT.draw_Gτ(G0_test, 0.1) == v[11]
    @test real(jDMFT.draw_Gτ(G0_test, 0.105)) > v[11] && real(jDMFT.draw_Gτ(G0_test, 0.105)) < v[12]
    @test jDMFT.draw_Gτ(G0_test, 2π) ≈ 0
    @test jDMFT.draw_Gτ(G0_test, -1.0) ≈ - jDMFT.draw_Gτ(G0_test, -1.0+2π)
    @test jDMFT.draw_Gτ(G0_test, -1.0 - 2*2π) ≈ jDMFT.draw_Gτ(G0_test, -1.0+2*2π)
end

@testset "tail subtraction" begin
    @test all(jDMFT.subtract_tail(1.0 ./ νnGrid, [0.0,1.0], νnGrid) .≈ 0.0)
end

@testset "Dyson Eq" begin
     GW, GImp, Σ = solve_AtomicLimit(νnGrid, μ, n, β, U)
     GLoc = jDMFT.GLoc(νnGrid, μ, kG, Σ)
     Δ    = jDMFT.Δ_FromGLocΣ(νnGrid, μ, GLoc, Σ)    
     Δc1  = jDMFT.kintegrate(kG, jDMFT.dispersion(kG) .^ 2)
     @test (- imag(GLoc) .* imag(νnGrid))[end] ≈ 1.0 atol=0.001
     @test (- imag(Δ) .* imag(νnGrid))[end] ≈ Δc1 atol=0.001

end
