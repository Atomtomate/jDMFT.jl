module jDMFT

using Random, LinearAlgebra, Logging
using FastGaussQuadrature
using FFTW
#using Dispersions
#TODO: tmporary workaround
const KGrid = Vector{Float64}

export sample
export MatsubaraFunction, τFunction
export subtract_tail
export τIntegrate, τ_to_ω, τ_to_ω, ω_to_τ

# impurity solvers
export impSolve_IPT

const DEBUG = false
struct NotImplemented <: Exception end

include("GFTools.jl")
include("DataTypes.jl")
include("Random.jl")
include("DMFTLoop.jl")
include("Fourier.jl")

include("IPT.jl")
include("CTINT.jl")

end # module
