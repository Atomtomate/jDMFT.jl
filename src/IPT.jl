# ==================================================================================================== #
#                                               IPT.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Julian Stobbe                                                                    #
#   Last Edit Date  : 20.08.22                                                                         #
# ----------------------------------------- Description ---------------------------------------------- #
#   Iterated perturbation theory impurity solver                                                       #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #

function impSolve_IPT(GWeiss_up::MatsubaraFunction, GWeiss_do::MatsubaraFunction, νnGrid::Vector{ComplexF64}, 
                      U::Float64, n::Float64)
    ΣImp_i = Vector{ComplexF64}(undef, length(GWeiss_up.data))
    impSolve_IPT!(ΣImp_i, GWeiss_up, GWeiss_do, νnGrid, U, n)
    return ΣImp_i
end

"""
    impSolve_IPT(GWeiss_up::τFunction, GWeiss_do::τFunction, νnGrid::Vector{ComplexF64}, U::Float64, n::Float64)

TODO: Documentation
"""
function impSolve_IPT!(Σ_new::Vector{ComplexF64}, GWeiss_up::MatsubaraFunction, GWeiss_do::MatsubaraFunction, νnGrid::Vector{ComplexF64}, 
                      U::Float64, n::Float64)
    Nτ = length(GWeiss_up.data)
    τGrid, τWeights = riemann(0.0,GWeiss_up.β-1/Nτ,Nτ)
    GW_τ_up = ω_to_τ(GWeiss_up, τWeights, τGrid)
    Σ_new[:] = ((U^2)/4) .* τ_to_ω(GW_τ_up.data .^ 3, νnGrid, GW_τ_up.τGrid, [0.0, 1.0], [0.5^3, 0.0], GW_τ_up.β)
end

function impSolve_IPT_MF(GWeiss_up::MatsubaraFunction, GWeiss_do::MatsubaraFunction, νnGrid::Vector{ComplexF64}, 
                      U::Float64, n::Float64)
    Σ = ((U^2)/4) .* τ_to_ω(GW_τ_up.data .^ 3, νnGrid, GW_τ_up.τGrid, [0.0, 1.0], [0.5^3, 0.0], GW_τ_up.β)
end
