# ================================================================================
#                   MatsubaraFunction and τFunction overload                     =
# ================================================================================

function  τ_to_ω(F::τFunction, νGrid::Vector{ComplexF64}, tail_coeffs::AbstractVector, ft_tail::AbstractVector)
    vals = τ_to_ω_GR(F.data, νGrid, F.τWeights, F.τGrid, F.β)
    MatsubaraFunction(vals, F.β, νGrid, tail_coeffs, ft_tail)
end

function  ω_to_τ(F::MatsubaraFunction, τWeights::AbstractVector, τGrid::AbstractVector)
    vals = ω_to_τ(F.data, F.νGrid, τGrid, F.tail_coeffs, F.ft_tail, F.β)
    τFunction(vals, F.β, τGrid, τWeights) 
end


# ================================================================================
#                           Low Level Implementations                            =
# ================================================================================

function τ_to_ω_GR(in::Vector{ComplexF64}, νnGrid::Vector{ComplexF64}, τWeights::Vector{Float64}, τGrid::Vector{Float64}, β::Float64)
    res = zeros(ComplexF64,length(νnGrid))
    w = (last(τGrid) - first(τGrid))/2
    for i in 1:length(res)
        νn = νnGrid[i]
        res[i] = β * sum(τWeights .* exp.(νn .* τGrid) .* in) * w
    end
    return res
end

function τ_to_ω(in::Vector{ComplexF64}, νnGrid::Vector{ComplexF64}, τGrid::Vector{Float64}, β::Float64)
    res = [β * sum(exp.(νn .* τGrid) .* in) / length(τGrid) for νn in νnGrid]
    return res
end

function ω_to_τ(in::Vector{ComplexF64}, νnGrid::Vector{ComplexF64}, τGrid::Vector{Float64},
                      tail_c::Vector{Float64}, tail_c_f::Vector{Float64}, β::Float64)
    res = zeros(ComplexF64,length(τGrid))
    in_sub = deepcopy(in)
    for i in 1:length(tail_c)
        subtract_tail!(in_sub, in_sub, tail_c[i], νnGrid .^ i)
    end

    for i in 1:length(res)
        for j in 1:length(νnGrid)
            νn = νnGrid[j]
            res[i] += exp(-νn * τGrid[i]) * in_sub[j]
        end
        res[i] = res[i]/β + sum(tail_c_f) 
    end
    return res
end

τIntegrate(f::Function, τWeights::Vector{Float64}, τGrid::Vector{Float64}) = dot(τWeights, f.(τGrid)) * (last(τGrid) - first(τGrid))/2
τIntegrate(f::τFunction) = dot(f.τWeights, f.data) * (last(f.τGrid) - first(f.τGrid))/2
