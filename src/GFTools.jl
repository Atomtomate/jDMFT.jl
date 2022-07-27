# ================================================================================
#                   Matsubara and imaginary time function-types                  =
# ================================================================================

struct MatsubaraFunction
    data::Vector{ComplexF64}
    β::Float64
    νGrid::Vector{ComplexF64}
    tail_coeffs::AbstractVector
    ft_tail::AbstractVector
end

struct τFunction
    data::Vector{ComplexF64}
    β::Float64
    τGrid::AbstractVector
    τWeights::AbstractVector
end


#TODO: test and implement fast-version (only draw random indices)
#TODO: replace linear interpolation by some cubic lin. inter, using legndre grid 
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

#TODO: no guaranteed to find nearest solution at index 1 and N
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

function subtract_tail!(outp::AbstractArray{T,1}, inp::AbstractArray{T,1}, c::Float64, iω::Array{ComplexF64,1}) where T <: Number
    for n in 1:length(inp)
        if iω[n] != 0
            outp[n] = inp[n] - c/iω[n]
        else
            outp[n] = inp[n]
        end
    end
end

function subtract_tail(inp::AbstractArray{T,1}, c::Float64, iω::Array{ComplexF64,1}) where T <: Number
    res = Array{eltype(inp),1}(undef, length(inp))
    subtract_tail!(res, inp, c, iω)
    return res
end


# ================================================================================
#                        Dyson Eq. related transformations                       =
# ================================================================================

Σ_Dyson(GBath::Array{ComplexF64,1}, GImp::Array{ComplexF64,1})::Array{ComplexF64,1} = 1 ./ GBath .- 1 ./ GImp

function GLoc(iνn::Vector{ComplexF64}, μ::Float64, kG::KGrid, ΣImp::Vector{ComplexF64})::Vector{ComplexF64}
    return [kintegrate(kG, 1 ./ (iν .+ μ .- kG.ϵkGrid .- ΣImp)) for iν in iνn]
end

WeissGF(GLoc::Vector{ComplexF64}, ΣImp::Vector{ComplexF64})::Vector{ComplexF64} = 1 ./ (1 ./ GLoc .+ ΣImp)

PH_transform(GWeiss::Vector{ComplexF64}, U::Float64) = 1 ./ (1 ./ GWeiss .- 0.5 * U)


# ================================================================================
#                                     Misc                                       =
# ================================================================================

iν_array(β::Real, grid::AbstractArray{Int64,1}) = ComplexF64[1.0im*((2.0 *el + 1)* π/β) for el in grid]
iν_array(β::Real, size::Integer)    = ComplexF64[1.0im*((2.0 *i + 1)* π/β) for i in 0:size-1]
iω_array(β::Real, grid::AbstractArray{Int64,1}) = ComplexF64[1.0im*((2.0 *el)* π/β) for el in grid]
iω_array(β::Real, size::Integer)    = ComplexF64[1.0im*((2.0 *i)* π/β) for i in 0:size-1]
