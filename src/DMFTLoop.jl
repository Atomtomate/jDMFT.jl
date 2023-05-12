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
    ΣImp_list = []
    tmp_Gup = deepcopy(GW_up)
    tmp_Gdo = deepcopy(GW_do)
    ΣImp_i = Vector{ComplexF64}(undef, length(GW_up.data))

    for i in 1:NIt

        ΣImp_i[:] = impSolve_IPT(tmp_Gup, tmp_Gdo, νnGrid, U, n)
        #push!(ΣImp_list, ΣImp_i)

        GLoc_i = GLoc(νnGrid, μ, kG, ΣImp_i)

        tmp_Gup = MatsubaraFunction(WeissGF(GLoc_i, ΣImp_i), β, νnGrid)
        tmp_Gdo = MatsubaraFunction(WeissGF(GLoc_i, ΣImp_i), β, νnGrid)
    end
    return tmp_Gup, ΣImp_i
end
