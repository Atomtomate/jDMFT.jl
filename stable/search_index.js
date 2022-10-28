var documenterSearchIndex = {"docs":
[{"location":"#jDMFT.jl-Documentation","page":"Home","title":"jDMFT.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Small demonstration of a Dynamical Mean Field Theory (DMFT) solver with some examples of impurity solvers included. For now there is a IteratedPerturbation Theory (IPT) and continuous time quantum Monte Carlo, interaction expansion (CT-INT), solver available.","category":"page"},{"location":"","page":"Home","title":"Home","text":"TODO: at some point I should provide some actual documentation here.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Autodocs","page":"Home","title":"Autodocs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [jDMFT]\nOrder   = [:module, :constant, :type, :function, :macro]","category":"page"},{"location":"#jDMFT.CTInt_Confs","page":"Home","title":"jDMFT.CTInt_Confs","text":"CTInt_Confs\n\nFields\n\nβ         : Float64, inverse temperature\nU         : Float64, interaction strenght\nδ         : Float64, numerical stabilization parameter for sample matrix\nGWeiss    : τFunction, noninteracting impurity Green's function\nτList     : Array{Float64}, τ values on which the Green's function is evaluated\nτiList    : Array{Float64}, sampled τ, i-th entry corresponds to i-th row/column in SampleMatrix\nsiList    : Array{Int}, external ising field, i-th entry corresponds to i-th row/column in SampleMatrix\n\n\n\n\n\n","category":"type"},{"location":"#jDMFT.MatsubaraFunction","page":"Home","title":"jDMFT.MatsubaraFunction","text":"MatsubaraFunction\n\nHolds data for a function over Matsubara frequencies and associated information. \n\nFields\n\ndata        : Vector{ComplexF64}, data\nβ           : Float64, inverse temperature\nfGrid       : Vector{ComplexF64}, x-values (i.e., frequence grid), frac(2n+1) pibeta for fermionic grids and frac(2n) pibeta for bosonic grids. \ntail_coeffs : AbstractVector, first N tail coefficients of the high frequency expansion sum_i=0^Nfracc_iinu^i_n.\n\n\n\n\n\n","category":"type"},{"location":"#jDMFT.Measurements","page":"Home","title":"jDMFT.Measurements","text":"Measurements\n\nAccumulation helper, see Gull et al., 2008.\n\n#TODO: this should be generalized to allow histograms and bootstrapping, e.g. dict of functors\n\nFields\n\nsamples      : Vector{ComplexF64}, List of all samples\nτGrid        : Vector{Float64}, sampling points for τ\nNSamples     : Int, number of samples\ntotalSign    : Int, total Monte Carlo sign.\ntotalExpOrder  : Int, sum of expansion orders\n\n\n\n\n\n","category":"type"},{"location":"#jDMFT.Measurements-Tuple{Int64, Float64, Symbol}","page":"Home","title":"jDMFT.Measurements","text":"Measurements(NBins::Int, β::Float64, type::Symbol)\n\nGenerates a Measurements struct for a temperature β with NBins bins and of grid type type which may be:\n\n'GaussQuad' : generates a Gauss Radau grid and its corresponding τWeights\n'Riemann'   : generates an evenly distriubted grid where each point has equal weight\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.SampleMatrix","page":"Home","title":"jDMFT.SampleMatrix","text":"SampleMatrix\n\nHolds samples data for f(τij) in each matrix entry. f depends on the algrithm (Weiss GF for CT-INT, Δ for CTHYB). This should be used in conjunction with the config struct (wich saves the τ_ij and additional data) for the specific sampler.\n\nFields\n\ndata        : Matrix{ComplexF64}, Data for sample matrix\nrowCache    : Vector{ComplexF64}, Row cache for potential new size\ncolCache    : Vector{ComplexF64}, Column cache for potential new size\nS           : Complex64, Cache for bottom right entry (remaining new value besides row and colum cache)\nN           : Int, Used matrix dimension, this differs from the actual data size, due to caching, see also update_cache!\n\n\n\n\n\n","category":"type"},{"location":"#jDMFT.τFunction","page":"Home","title":"jDMFT.τFunction","text":"τFunction\n\nHolds data for a function over imaginary time and associated information. \n\nFields\n\ndata        : Vector{ComplexF64}, data\nβ           : Float64, inverse temperature\nτGrid       : AbstractVector{Float64}, tau-grid, \nτWeights    : AbstractVector{Float64}, wheights for each f(tau) point, used for integration methods.\ntail_coeffs : AbstractVector, first N tail coefficients of the high frequency expansion sum_i=0^Nfracc_itau.\n\n\n\n\n\n","category":"type"},{"location":"#jDMFT.accumulate!-Tuple{Random.AbstractRNG, jDMFT.Measurements, Int64, jDMFT.CTInt_Confs, jDMFT.SampleMatrix}","page":"Home","title":"jDMFT.accumulate!","text":"accumulate!(rng::AbstractRNG, m::Measurements, sign::Int, confs::CTInt_Confs, M::SampleMatrix; with_τshift=true)\n\nUpdates Monte Carlo estimate by measuring impurity Green's function for given configuration stored in M.\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.draw_Gτ-Tuple{τFunction, Float64}","page":"Home","title":"jDMFT.draw_Gτ","text":"draw_Gτ(Gτ::τFunction, τi::Float64)::ComplexF64 \ndraw_Gτ(data::Vector{ComplexF64}, τi::Float64, τGrid::Vector{Float64}, β::Float64)::ComplexF64\n\nObtain f(\\tau) data for a function, either given as τFunction or as Vector{ComplexF64} at  imaginary time τi. Data is linearly interpolated from the given grid.\n\nTODO: test and implement fast-version (only draw random indices) TODO: replace linear interpolation by some cubic lin. inter, using legndre grid \n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.fast_update_decr!-Tuple{jDMFT.SampleMatrix, Int64}","page":"Home","title":"jDMFT.fast_update_decr!","text":"fast_update_decr!(M::SampleMatrix, N::Int)\n\nCompute fast inverse update, using Woodbury matrix idenity:      M_inv = inv(M)     Mnew_inv =M_inv after deleting right column and bottom row     Mnew = inv(Mnew_inv)\n\nN is the new size. TODO: currently only rank 1 updates supported, so N = M.N-1.\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.fast_update_incr!-Tuple{jDMFT.SampleMatrix, Int64}","page":"Home","title":"jDMFT.fast_update_incr!","text":"fast_update_incr!(M::SampleMatrix, N::Int)\n\nCompute fast inverse update, using Woodbury matrix idenity:      M_inv = inv(M)     Mnew_inv = M_inv after inserting row/col/corner obtained from M.rowCach, M.colCache, M.S     Mnew = inv(Mnew_inv)\n\nM.rowCache, M.colCache and M.S must be set to the new values, before calling this function! N is the new size. TODO: currently only rank 1 updates supported, so N = M.N+1.\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.gaussradau-Tuple{Real, Real, Int64}","page":"Home","title":"jDMFT.gaussradau","text":"gaussradau(start::Real, stop::Real, N::Int)\n\nGenerate grid for gaussradau sum integration from start to stop with N points. Returns (weights, grid), this is a wrapper around gaussradau from FastGaussQuadrature.\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.impSolve_IPT-Tuple{τFunction, τFunction, Vector{ComplexF64}, Float64, Float64}","page":"Home","title":"jDMFT.impSolve_IPT","text":"impSolve_IPT(GWeiss_up::τFunction, GWeiss_do::τFunction, νnGrid::Vector{ComplexF64}, U::Float64, n::Float64)\n\nTODO: Documentation\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.measure_GImp_τ-Tuple{jDMFT.Measurements, τFunction}","page":"Home","title":"jDMFT.measure_GImp_τ","text":"measure_GImp_τ(m::Measurements, GWeiss::τFunction)\n\nMeasure impurity Green's function, after Monte Carlo samples have been accumulated in m, see also Measurements.\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.riemann-Tuple{Real, Real, Int64}","page":"Home","title":"jDMFT.riemann","text":"riemann(start::Real, stop::Real, N::Int)\n\nGenerate grid for Riemann sum integration from start to stop with N points. Returns (weights, grid)\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.subtract_tail!-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}, Vector{Float64}, Vector{ComplexF64}}} where T<:Number","page":"Home","title":"jDMFT.subtract_tail!","text":"subtract_tail(inp::AbstractVector{T}, c::Vector{Float64}, iω::Vector{ComplexF64})\nsubtract_tail!(outp::AbstractVector{Number}, inp::AbstractVector{Number}, c::Vector{Float64}, iω::Vector{ComplexF64})\n\nSubtract high frequency tail of function, i.e. f(iomega_n) - sum_l fracc_liomega_n^l, with tail coefficients c_l and tail iω. One can use iν_array or iω_array to generate the grid. the inplace version stores the resulting data in outp.\n\nTODO: this function is not optimized for performance\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.τGridTransform-Tuple{Real, Real, Vector{Float64}}","page":"Home","title":"jDMFT.τGridTransform","text":"τGridTransform(start::Float64, stop::Float64, grid::Vector{Float64})\n\nTransform integration grid from [-1,1] to [start, stop].\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.τIndex-Tuple{Float64, Vector{Float64}, Float64}","page":"Home","title":"jDMFT.τIndex","text":"τIndex(τ::Float64, τGrid::Vector{Float64}, β::Float64)::Tuple{Int,Int}\n\nInternal function to determine nearest imaginary time value for τGrid of Green's function.  Find index for given τ that matches nearest value in τGrid. Returns Tuple with sign (since ``G(\\tau + \\beta) = - G(\\tau)) and index. \n\nTODO: no guaranteed to find nearest solution at index 1 and N\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.τIntegrate-Tuple{AbstractVector, Vector, Vector{Float64}}","page":"Home","title":"jDMFT.τIntegrate","text":"τIntegrate(f::AbstractVector, τWeights::Vector, τGrid::Vector{Float64})\nτIntegrate(f::Function, τWeights::Vector, τGrid::Vector{Float64})\nτIntegrate(f::τFunction)\n\nIntegrates Vector of values, Function or τFunction over τGrid, given τWeights.  Note, that you hav to transform any τ-spacing from the [-1,1] interval using τGridTransform.\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.τ_to_ω-Tuple{τFunction, Vector{ComplexF64}}","page":"Home","title":"jDMFT.τ_to_ω","text":"τ_to_ω(F::τFunction, fGrid::Vector{ComplexF64})\n\nFourier transform from imaginary time to Matsubara space with a grid given bei fGrid.\n\n\n\n\n\n","category":"method"},{"location":"#jDMFT.ω_to_τ-Tuple{MatsubaraFunction, AbstractVector, AbstractVector}","page":"Home","title":"jDMFT.ω_to_τ","text":"ω_to_τ(F::MatsubaraFunction, τWeights::AbstractVector, τGrid::AbstractVector)\n\nFourier transform from Matsubara frequencies to imaginary time with given τGrid and τWeights (for integration such as Gauss quadrature. See τGridTransform, gaussradau and riemann).\n\n\n\n\n\n","category":"method"}]
}