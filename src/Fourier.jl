# ==================================================================================================== #
#                                           Fourier.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Julian Stobbe                                                                    #
#   Last Edit Date  : 27.10.22                                                                         #
# ----------------------------------------- Description ---------------------------------------------- #
#   Fourier transform and integration routines.                                                        #
# -------------------------------------------- TODO -------------------------------------------------- #
#   - implement low level FFT capabilities.                                                            #
#   - implement low level MatsubaraSum function, similar to τIntegrate                                 #
#   - replace @asserts by throws statements                                                            #
# ==================================================================================================== #


# ================================================================================
# =                 MatsubaraFunction and τFunction overload                     =
# ================================================================================

"""
    τ_to_ω(F::τFunction, fGrid::Vector{ComplexF64})

Fourier transform from imaginary time to Matsubara space with a grid given bei `fGrid`.
"""
function τ_to_ω(F::τFunction, fGrid::Vector{ComplexF64})
    ω_coeffs = tailTransform_τ_to_ω(F.tail_coeffs)
    vals = τ_to_ω(F.data, fGrid, F.τGrid, ω_coeffs, F.tail_coeffs, F.β)
    MatsubaraFunction(vals, F.β, fGrid, tailTransform_τ_to_ω(F.tail_coeffs))
end

"""
    ω_to_τ(F::MatsubaraFunction, τWeights::AbstractVector, τGrid::AbstractVector)

Fourier transform from Matsubara frequencies to imaginary time with given `τGrid` and `τWeights` (for integration such as Gauss quadrature. See [`τGridTransform`](@ref), [`gaussradau`](@ref) and [`riemann`](@ref)).
"""
function ω_to_τ(F::MatsubaraFunction, τWeights::AbstractVector, τGrid::AbstractVector)
    τ_coeffs = tailTransform_ω_to_τ(F.tail_coeffs)
    vals = ω_to_τ(F.data, F.fGrid, τGrid, F.tail_coeffs, τ_coeffs, F.β)
    τFunction(vals, F.β, τGrid, τWeights, τ_coeffs) 
end


# ================================================================================
# =                         Low Level Implementations                            =
# ================================================================================

function τ_to_ω(data_τ::Vector{ComplexF64}, νnGrid::Vector{ComplexF64}, τGrid::Vector{Float64},
                  ω_coeffs::Vector{Float64}, τ_coeffs::Vector{Float64}, β::Float64)
    @assert length(data_τ) == length(τGrid)
    prep_arr  = β .* exp.(1im * π * (k/length(τGrid)) for k in 0:(length(τGrid)-1)) .* (data_τ .+ τ_coeffs[1])
    return subtract_tail(fftshift(ifft(prep_arr)), -1.0 .* ω_coeffs, νnGrid)
end

function ω_to_τ(data_ω::Vector{ComplexF64}, νnGrid::Vector{ComplexF64}, τGrid::Vector{Float64},
                      ω_coeffs::Vector{Float64}, τ_coeffs::Vector{Float64}, β::Float64)
    @assert length(data_ω) == length(νnGrid)

    Nν  = (imag(first(νnGrid))*β/π - 1) / 2 
    inp = subtract_tail(data_ω, ω_coeffs,νnGrid)
    tmp = fft(inp) ./ β
    return tmp .* exp.(-2 * π * 1im * (-Nν + 1/2) .* [k/length(τGrid) for k in 0:length(τGrid)-1] ) .- τ_coeffs[1]
end


# ================================================================================
# =                               Grid Functions                                 =
# ================================================================================

"""
    riemann(start::Real, stop::Real, N::Int)

Generate grid for Riemann sum integration from `start` to `stop` with `N` points.
Returns `(weights, grid)`
"""
riemann(start::Real, stop::Real, N::Int) = τGridTransform(start, stop, collect(LinRange(-1, 1, N))), 2 .*ones(N) ./ N

"""
    gaussradau(start::Real, stop::Real, N::Int)

Generate grid for gaussradau sum integration from `start` to `stop` with `N` points.
Returns `(weights, grid)`, this is a wrapper around `gaussradau` from `FastGaussQuadrature`.
"""
function gaussradau(start::Real, stop::Real, N::Int) 
    grid,weights = FastGaussQuadrature.gaussradau(N)
    τGridTransform(start, stop, grid), weights
end

"""
    τGridTransform(start::Float64, stop::Float64, grid::Vector{Float64})

Transform integration grid from `[-1,1]` to `[start, stop]`.
"""
τGridTransform(start::Real, stop::Real, grid::Vector{Float64}) = ((stop - start) / 2) .* grid .+ (start + stop) / 2

"""
    τIntegrate(f::AbstractVector, τWeights::Vector, τGrid::Vector{Float64})
    τIntegrate(f::Function, τWeights::Vector, τGrid::Vector{Float64})
    τIntegrate(f::τFunction)

Integrates `Vector` of values, `Function` or `τFunction` over `τGrid`, given `τWeights`. 
Note, that you hav to transform any τ-spacing from the `[-1,1]` interval using [`τGridTransform`](@ref).
"""
τIntegrate(f::AbstractVector, τWeights::Vector, τGrid::Vector{Float64}) = dot(τWeights, f) * (last(τGrid) - first(τGrid))/2
τIntegrate(f::Function, τWeights::Vector, τGrid::Vector{Float64}) = dot(τWeights, f.(τGrid)) * (last(τGrid) - first(τGrid))/2
τIntegrate(f::τFunction) = dot(f.τWeights, f.data) * (last(f.τGrid) - first(f.τGrid))/2

