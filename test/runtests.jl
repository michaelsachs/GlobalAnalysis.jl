using Test

@testset "GlobalAnalysis.jl" begin
    include("jupyterNotebooks.jl")
    include("kineticFit.jl")
end
