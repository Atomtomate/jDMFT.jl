function init_weissgf_from_ed(fn::String, Niν::Int)
    N, μ, β, ϵp, tp  = 0, 0, 0, Float64[], Float64[]
    i = 0
    lines = readlines(fn)
    for l in lines
        if i == 4
            bb = split(l)[1]
            β = parse(Float64, bb[1:findfirst("d", bb)[1]-1])
        elseif i == 6
            N  = parse(Int, split(l,",")[1])
        elseif i >= 9 && i < 9+N-1
            push!(ϵp, parse(Float64, l))
        elseif i > 9+N-1 && i < 9+2*N-1
            push!(tp, parse(Float64, l))
        elseif i >= 9+2*N-3
            μ = parse(Float64, split(l)[1])
        end
        i += 1
    end
    #println("Read from file: N=", N, ", μ=", μ, ", β=", β, ", ϵp=", ϵp, ", tp=", tp)
    init_weissgf_from_ed(Niν, μ, β, ϵp, tp)
end

function init_weissgf_from_ed(N::Int, μ::Float64, β::Float64, ϵp::Vector{Float64}, tp::Vector{Float64})
    iν_n = jDMFT.iν_array(β, -N:N-1)
    Δ = [ sum( (tp .* conj(tp)) ./ (iν .- ϵp)) for iν in iν_n]
    G = 1 ./ (iν_n .+ μ .- Δ )
    return G
end

function init_weissgf_U0guess(N::Int, μ::Float64, β::Float64)
    ϵp = [0.5, 1.0, -0.5, -1.0]
    tp = [0.25, 0.25, 0.25, 0.25]
    init_weissgf_from_ed(N, μ, β, ϵp, tp)
end

function init_GImp(fn::String, Nν::Int)
    arr = readdlm(fn)
    νn = arr[1:Nν,1]
    β  = π/νn[1]
    data = arr[1:Nν,2] .+ 1im .* arr[1:Nν,3]
    νn = [-1 .* reverse(νn[2:end]); νn]
    G = MatsubaraFunction([conj(reverse(data[2:end])); data], β, 1im .* νn, [1.0], [-0.5])
    return G
end

function solve_AtomicLimit(νnGrid::Vector{ComplexF64}, μ::Float64, n::Float64, β::Float64, U::Float64)
    βi = BigFloat(β)
    μi = BigFloat(μ)
    Ui = BigFloat(U)
    n = Float64(2 * (exp(βi*μi) + exp(βi*(2*μi-Ui))) / (1 + 2*exp(βi*μi) + exp(βi*(2*μi-Ui))))
    GW    = [1 / (iν + μ) for iν in νnGrid]
    GImp  = [(1 - n/2) / (iν + μ) + (n/2)/(iν + μ - U) for iν in νnGrid]
    Σ     = 1 ./ GW .- 1 ./ GImp
    return GW, GImp, Σ
end

function solve_NonIntLimit(νnGrid::Vector{ComplexF64}, kG::KGrid, μ::Float64, n::Float64, β::Float64, U::Float64)

    Σ    = zeros(ComplexF64, length(νnGrid))
    GLoc = jDMFT.GLoc(νnGrid, μ::Float64, kG, Σ)
    GImp = deepcopy(GLoc)
    GW   = deepcopy(GLoc)
    return GW, GImp, Σ
end
