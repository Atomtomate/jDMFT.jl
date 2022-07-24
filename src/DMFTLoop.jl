function DMFTLoop(G0, kG::KGrid, μ::Float64, β::Float64)
    GImp = ImpSolve(G0, β)
    ΣImp_i = Σ_Dyson(G0, GImp)
    GLoc_i = GLoc(iνn, μ, kG, ΣImp)
    G0 = WeissGF(GLoc_i, ΣImp_i)
    return ΣImp_i
end
