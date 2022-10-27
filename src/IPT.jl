# ==================================================================================================== #
#                                               IPT.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Julian Stobbe                                                                    #
#   Last Edit Date  : 20.08.22                                                                         #
# ----------------------------------------- Description ---------------------------------------------- #
#   Iterated perturbation theory impurity solver                                                       #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #


"""
    impSolve_IPT(GWeiss_up::τFunction, GWeiss_do::τFunction, νnGrid::Vector{ComplexF64}, U::Float64, n::Float64)

TODO: Documentation
"""
function impSolve_IPT(GWeiss_up::τFunction, GWeiss_do::τFunction, νnGrid::Vector{ComplexF64}, 
                      U::Float64, n::Float64)
    Σ = U*n/2 .+ U^2 .* τ_to_ω(GWeiss_up.data .^ 3, νnGrid, GWeiss_up.τGrid, [0.5^3, 0.0], [0.0, 0.5^2], GWeiss_up.β)
    return Σ   
end
