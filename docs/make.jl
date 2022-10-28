push!(LOAD_PATH,"../src/")
using Documenter
using Pkg
Pkg.activate(String(@__DIR__) * "/..")
using jDMFT

makedocs(sitename="jDMFT")

makedocs(;
    modules=[jDMFT],
    authors="Julian Stobbe <Atomtomate@gmx.de> and contributors",
    repo="https://github.com/Atomtomate/jDMFT.jl/blob/{commit}{path}#L{line}",
    sitename="jDMFT",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://Atomtomate.github.io/jDMFT.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    branch="gh-pages",
    devbranch = "master",
    devurl = "stable",
    repo="github.com/Atomtomate/jDMFT.jl.git",
)
