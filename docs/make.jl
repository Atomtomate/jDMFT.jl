push!(LOAD_PATH,"../src/")
using Documenter
using Pkg
Pkg.activate(String(@__DIR__) * "/..")
using jDMFT

makedocs(sitename="jDMFT")
