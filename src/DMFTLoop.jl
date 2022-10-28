# ==================================================================================================== #
#                                           DMFTLoop.jl                                                #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Julian Stobbe                                                                    #
#   Last Edit Date  : 28.10.22                                                                         #
# ----------------------------------------- Description ---------------------------------------------- #
#   DMFT loop implementation                                                                           #
# -------------------------------------------- TODO -------------------------------------------------- #
#   - This is only a stub for now.                                                                     #
# ==================================================================================================== #



function DMFTLoop(NIt::Int, GW_up::MatsubaraFunction, GW_do::MatsubaraFunction, kG::KGrid, U::Float64, n::Float64, μ::Float64, β::Float64)
    νnGrid = GW_up.fGrid
    Nτ     = length(νnGrid)
    τGrid, τWeights = riemann(0.0,β-1/Nτ,Nτ)
    ΣImp_list = []

    for i in 1:NIt

        GW_τ_up = ω_to_τ(GW_up, τWeights, τGrid)
        GW_τ_do = ω_to_τ(GW_do, τWeights, τGrid)

        ΣImp_i = impSolve_IPT(GW_τ_up, GW_τ_do, νnGrid, U, n)
        push!(ΣImp_list, ΣImp_i)

        GLoc_i = GLoc(νnGrid, μ, kG, ΣImp_i)

        GW_up = MatsubaraFunction(WeissGF(GLoc_i, ΣImp_i), β, νnGrid)
        GW_do = MatsubaraFunction(WeissGF(GLoc_i, ΣImp_i), β, νnGrid)
    end
    return ΣImp_list
end
