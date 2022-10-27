# ==================================================================================================== #
#                                           GFTools.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Julian Stobbe                                                                    #
#   Last Edit Date  : 26.10.22                                                                         #
# ----------------------------------------- Description ---------------------------------------------- #
#   Green's function tools and data types.                                                             #
# -------------------------------------------- TODO -------------------------------------------------- #
#   - This should be a independent module (exists as file in 4+ projects).                             #
#   - several improvements to functions and separation into logical units is missing                   #
#   - Tail coefficients should be made consistent                                                      #
#       all types should know their 0,1/x,1/x^2.. tails                                                #
#       independent (global) lookup for analyticsums/FT                                                #
#   - Implement bosonic grid types                                                                     #
# ==================================================================================================== #


# ================================================================================
#                               Tail transform helpers                           =
# ================================================================================

function tailTransform_ω_to_τ(tail_coeffs::AbstractArray)
    length(tail_coeffs) != 2 && throw("Only length 2 tail coefficient transformations implemented yet.")
    tail_coeffs[1] != 0      && throw("constant coefficient of Matsubara tail needs to be 0!")
    return [0.5*tail_coeffs[2], 0.0]
end

function tailTransform_τ_to_ω(tail_coeffs::AbstractArray) 
    length(tail_coeffs) != 2 && throw("Only length 2 tail coefficient transformations implemented yet.")
    tail_coeffs[2] != 0      && throw("1/τ coefficient of imaginary time tail needs to be 0!")
    return [0.0, 2.0*tail_coeffs[1]]
end


# ================================================================================
#                   Matsubara and imaginary time function-types                  =
# ================================================================================

"""
    MatsubaraFunction

Holds data for a function over Matsubara frequencies and associated information. 

Fields
-------------
- **`data`**        : `Vector{ComplexF64}`, data
- **`β`**           : `Float64`, inverse temperature
- **`fGrid`**       : `Vector{ComplexF64}`, x-values (i.e., frequence grid), ``\\frac{(2n+1) \\pi}{\\beta}`` for fermionic grids and ``\\frac{(2n) \\pi}{\\beta}`` for bosonic grids. 
- **`tail_coeffs`** : `AbstractVector`, first `N` tail coefficients of the high frequency expansion ``\\sum_{i=0}^N\\frac{c_i}{i\\nu^i_n}``.
"""
struct MatsubaraFunction
    data::Vector{ComplexF64}
    β::Float64
    fGrid::Vector{ComplexF64}
    tail_coeffs::AbstractVector
end

"""
    τFunction

Holds data for a function over imaginary time and associated information. 

Fields
-------------
- **`data`**        : `Vector{ComplexF64}`, data
- **`β`**           : `Float64`, inverse temperature
- **`τGrid`**       : `AbstractVector{Float64}`, ``\\tau``-grid, 
- **`τWeights`**    : `AbstractVector{Float64}`, wheights for each ``f(\\tau)`` point, used for integration methods.
- **`tail_coeffs`** : `AbstractVector`, first `N` tail coefficients of the high frequency expansion ``\\sum_{i=0}^N\\frac{c_i}{\\tau}``.
"""
struct τFunction
    data::Vector{ComplexF64}
    β::Float64
    τGrid::AbstractVector{Float64}
    τWeights::AbstractVector{Float64}
    tail_coeffs::AbstractVector
end


"""
    draw_Gτ(Gτ::τFunction, τi::Float64)::ComplexF64 
    draw_Gτ(data::Vector{ComplexF64}, τi::Float64, τGrid::Vector{Float64}, β::Float64)::ComplexF64

Obtain `f(\\tau)` data for a function, either given as [`τFunction`](@ref) or as `Vector{ComplexF64}` at 
imaginary time `τi`. Data is linearly interpolated from the given grid.

TODO: test and implement fast-version (only draw random indices)
TODO: replace linear interpolation by some cubic lin. inter, using legndre grid 
"""
draw_Gτ(Gτ::τFunction, τi::Float64)::ComplexF64 = draw_Gτ(Gτ.data, τi, Gτ.τGrid, Gτ.β)

function draw_Gτ(data::Vector{ComplexF64}, τi::Float64, τGrid::Vector{Float64}, β::Float64)::ComplexF64
    Nτ  = length(data)
    n    = floor(Int, τi/β)
    sign = 1 - 2*(abs(n) % 2)
    τi    = τi - n*β
    ind  = searchsortedfirst(τGrid, τi)
    sign, ind
    ind = ind > Nτ ? ind-1 : ind
    m = if ind < Nτ 
        (data[ind+1] - data[ind])/(τGrid[ind+1] - τGrid[ind]) 
    else
        (data[1] - data[ind])/(τGrid[1] - τGrid[ind] + β)
    end
    val = data[ind] + m * (τi - τGrid[ind])
    return sign*val
end


"""
    τIndex(τ::Float64, τGrid::Vector{Float64}, β::Float64)::Tuple{Int,Int}

Internal function to determine nearest imaginary time value for `τGrid` of Green's function. 
Find index for given `τ` that matches nearest value in `τGrid`. Returns `Tuple` with
sign (since ``G(\\tau + \\beta) = - G(\\tau)) and index. 

TODO: no guaranteed to find nearest solution at index 1 and N
"""
function τIndex(τ::Float64, τGrid::Vector{Float64}, β::Float64)::Tuple{Int,Int}
    n    = floor(Int, τ/β)
    sign = 1 - 2*(abs(n) % 2)
    τ    = τ - n*β
    ind  = searchsortedfirst(τGrid, τ)
    ind  = ind > length(τGrid) ? ind - 1 : ind
    sign, ind
end


# ================================================================================
#                               Tail subtraction                                 =
# ================================================================================

"""
    subtract_tail(inp::AbstractVector{T}, c::Vector{Float64}, iω::Vector{ComplexF64})
    subtract_tail!(outp::AbstractVector{Number}, inp::AbstractVector{Number}, c::Vector{Float64}, iω::Vector{ComplexF64})

Subtract high frequency tail of function, i.e. ``f(i\\omega_n) - \\sum_l \\frac{c_l}{i\\omega_n^l}``, with tail
coefficients `c_l` and tail `iω`.
One can use [`iν_array`](@ref) or [`iω_array`](@ref) to generate the grid.
the inplace version stores the resulting data in `outp`.

TODO: this function is not optimized for performance
"""
function subtract_tail!(outp::AbstractVector{T}, inp::AbstractVector{T}, c::Vector{Float64}, iω::Vector{ComplexF64}) where T <: Number
    outp[:] = deepcopy(inp) .- c[1]
    for (l,cl) in enumerate(c[2:end])
        for n in 1:length(inp)
            if iω[n] != 0
                outp[n] = outp[n] - cl/(iω[n]^l)
            else
                outp[n] = outp[n]
            end
        end
    end
end

function subtract_tail(inp::AbstractVector{T}, c::Vector{Float64}, iω::Vector{ComplexF64}) where T <: Number
    res = Vector{eltype(inp)}(undef, length(inp))
    subtract_tail!(res, inp, c, iω)
    return res
end


# ================================================================================
#                        Dyson Eq. related transformations                       =
# ================================================================================

Σ_Dyson(GBath::Vector{ComplexF64}, GImp::Vector{ComplexF64})::Vector{ComplexF64} = 1 ./ GBath .- 1 ./ GImp

function GLoc(iνn::Vector{ComplexF64}, μ::Float64, kG::KGrid, ΣImp::Vector{ComplexF64})::Vector{ComplexF64}
    return [kintegrate(kG, 1 ./ (iν .+ μ .- kG.ϵkGrid .- ΣImp)) for iν in iνn]
end

WeissGF(GLoc::Vector{ComplexF64}, ΣImp::Vector{ComplexF64})::Vector{ComplexF64} = 1 ./ (1 ./ GLoc .+ ΣImp)

PH_transform(GWeiss::Vector{ComplexF64}, U::Float64, α::Float64=0.5) = 1 ./ (1 ./ GWeiss .- α * U)


# ================================================================================
#                                     Misc                                       =
# ================================================================================

iν_array(β::Real, grid::AbstractVector{Int64}) = ComplexF64[1.0im*((2.0 *el + 1)* π/β) for el in grid]
iν_array(β::Real, size::Integer)    = ComplexF64[1.0im*((2.0 *i + 1)* π/β) for i in 0:size-1]
iω_array(β::Real, grid::AbstractVector{Int64}) = ComplexF64[1.0im*((2.0 *el)* π/β) for el in grid]
iω_array(β::Real, size::Integer)    = ComplexF64[1.0im*((2.0 *i)* π/β) for i in 0:size-1]
