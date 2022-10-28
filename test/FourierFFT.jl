using FFTW

#νnGrid = jDMFT.iν_array(β, -Nν:Nν-1); 
data_τ = G_τ.data
grid_τ = G_τ.τGrid
data_ω = GW.data

# ==================================================================================================== #
#                                               ω to τ                                                 #
# ==================================================================================================== #

res_τ_sum = zeros(ComplexF64,length(grid_τ))
data_ω = GW.data
data_ω_sub = deepcopy(data_ω)
tail_c = [1.0]
tail_c_f = [0.5]
for i in 1:length(tail_c)
    for n in 1:length(data_ω)
        data_ω_sub[n] = data_ω[n] - tail_c[i]/νnGrid[n]
    end
end
for i in 1:length(res_τ_sum)
    for (j,νn) in enumerate(νnGrid)
        res_τ_sum[i] += exp(-νn * grid_τ[i]) * data_ω_sub[j]
    end
    res_τ_sum[i] = res_τ_sum[i]/β + sum(tail_c_f) 
end

inp = data_ω_sub
tmp = fft(inp) ./ β
res_τ_new = tmp .* exp.(-2 * π * 1im * (-Nν + 1/2) .* [k/length(grid_τ) for k in 0:length(grid_τ)-1] ) .- 0.5
 

# ==================================================================================================== #
#                                               τ to ω                                                 #
# ==================================================================================================== #

data_τ = deepcopy(res_τ_new)

res_sum = zeros(ComplexF64, length(νnGrid))
for (i,νn) in enumerate(νnGrid)
    for (j,τ) in enumerate(grid_τ)
        res_sum[i] += β/length(grid_τ) * exp(νn*τ) * data_τ[j]
    end
end


prep_arr  = β .* exp.(1im * π * (k/length(grid_τ)) for k in 0:(length(grid_τ)-1)) .* (data_τ .+ 0.5)
res_new   = (fftshift(ifft(prep_arr)) ) .+ 1 ./ νnGrid

# println("real parts check0 : ", round(sum(abs.(real(GW.data))), digits=6))
# println("           sum    : ", round(sum(abs.(real(res_sum))),digits=6))
# println("           new    : ", round(sum(abs.(real(res_new))),digits=6))
# println("imag dif   check0 : ", round(sum(imag(GW.data .- GW.data)), digits=6))
# println("           sum    : ", round(sum(abs.(imag(res_sum .- GW.data))),digits=6))
# println("           new    : ", round(sum(abs.(imag(res_new .- GW.data))),digits=6))

@testset "FFT Fourier" begin
    @test sum(abs.(real(GW.data))) ≈ 0 atol=1e-8
    @test sum(abs.(real(res_new))) ≈ 0 atol=1e-8
    @test sum(imag(GW.data .- GW.data)) ≈ 0 atol=1e-8
    @test sum(abs.(imag(res_new .- GW.data))) ≈ 0 atol=1e-8
end
