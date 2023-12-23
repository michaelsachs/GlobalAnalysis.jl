using Documenter
using GlobalAnalysis

makedocs(
    sitename = "GlobalAnalysis.jl",
    authors = "Michael Sachs and contributors.",
    format = Documenter.HTML(),
    modules = [GlobalAnalysis],

    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Data format" => "import.md",
        "Kinetic model" => "kineticModel.md",
    ],
)

deploydocs(
    repo = "github.com/michaelsachs/GlobalAnalysis.jl.git",
    devbranch = "main"
)