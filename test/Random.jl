@testset "Random Number Generator" begin
    @test typeof(jDMFT.rng) <: MersenneTwister
end
