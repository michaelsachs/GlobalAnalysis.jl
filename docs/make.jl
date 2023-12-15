using Pkg
pkg"activate .."
push!(LOAD_PATH,"../src/")

using Documenter
using GlobalAnalysis

makedocs(
    sitename = "GlobalAnalysis",
    authors = "Michael Sachs and contributors.",
    format = Documenter.HTML(),
    modules = [GlobalAnalysis],

    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Data import" => "import.md",
        "Kinetic model" => "kineticModel.md",
    ],
)

deploydocs(
    repo = "github.com/michaelsachs/GlobalAnalysis.jl.git",
)