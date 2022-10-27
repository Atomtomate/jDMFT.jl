@testset "SampleMatrix" begin
    M0 = jDMFT.SampleMatrix()
    M1 = jDMFT.SampleMatrix(15)
    @test typeof(M0.data) <: Matrix
    @test M1.N == 0
    @test all(size(M1.data) .== (15,15))
    @test all(size(M1.colCache) .== (15,))
    @test all(size(M1.rowCache) .== (15,))
end

@testset "fast_update_incr!" begin
    Mt = 1 ./ collect(reshape(1:2.0:18,3,3)) # invertible 3x3 matrix
    M_inv = inv(Mt)
    M = jDMFT.SampleMatrix(5)
    M.data[1:3,1:3] .= deepcopy(Mt)
    M.N              = 3                     # fake 3x3 state

    Mnew_inv = Matrix{Float64}(undef, 4,4)   # new invertex matrix
    Mnew_inv[1:3,1:3] .= M_inv
    Mnew_inv[4,1:4] = collect(11.1:1.4:16)
    Mnew_inv[1:3,4] = collect(5.1:1.4:8)
    Mnew  = inv(Mnew_inv)                    # New matrix

    M.rowCache = deepcopy(Mnew_inv[4,1:4])
    M.colCache = deepcopy(Mnew_inv[1:3,4])
    M.S = Mnew_inv[4,4] .+ 0.0im
    jDMFT.fast_update_incr!(M, 4)
    @test all(M.data[1:4,1:4] .≈ Mnew)
end

@testset "fast_update_decr!" begin
    Mt = 1 ./ collect(reshape(1:2.0:18,3,3)) # invertible 3x3 matrix
    M_inv = inv(Mt)
    M = jDMFT.SampleMatrix(5)
    M.data[1:3,1:3] .= deepcopy(Mt)
    M.N              = 3                     # fake 3x3 state

    Mnew  = inv(M_inv[1:2,1:2])

    jDMFT.fast_update_decr!(M, 2)
    @test all(M.data[1:2,1:2] .≈ Mnew)
end
