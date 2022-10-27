
function DMFTLoop(NIt::Int, GW_up::MatsubaraFunction, GW_τ_do::MatsubaraFunction, kG::KGrid, U::Float64, n::Float64, μ::Float64, β::Float64)
    D      = 1.0
    νnGrid = GW_up.fGrid
    Nτ     = length(νnGrid)
    τGrid, τWeights = riemann(0.0,β-1/Nτ,Nτ)
    ΣImp_list = []

    for i in 1:NIt

        GW_τ_up = ω_to_τ(GW_up, τWeights, τGrid)
        GW_τ_do = ω_to_τ(GW_do, τWeights, τGrid)

        ΣImp_i = impSolve_IPT(GW_τ_up, GW_τ_do, νnGrid, U, n)
        push!(ΣImp_list, ΣImp_i)

        # ΣImp_i_up = Σ_Dyson(GW_up, GImp)
        # ΣImp_i_do = Σ_Dyson(GW_do, GImp)
        #
        #GLoc_i = GLoc(iνn, μ, kG, ΣImp)
        #G0 = WeissGF(GLoc_i, ΣImp_i)

        GImp = GW_up.data ./ (1 .- GW_up.data .* ΣImp_i)
        GW_up = MatsubaraFunction(1 ./ (νnGrid .- D .* D .* GImp ./ 4), GW_up.β, νnGrid, [0.0, 1.0])
        GW_do = MatsubaraFunction(1 ./ (νnGrid .- D .* D .* GImp ./ 4), GW_do.β, νnGrid, [0.0, 1.0])
    end
    return ΣImp_list
end
