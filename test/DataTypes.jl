@testset "SampleMatrix" begin
    M0 = SampleMatrix()
    M1 = SampleMatrix(15)
    @test typeof(M0.data) <: Matrix
    @test M1.N == 0
    @test all(size(M1.data) .== (15,15))
    @test all(size(M1.colCache) .== (15,))
    @test all(size(M1.rowCache) .== (15,))
end
