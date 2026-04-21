# set BLAS threads to 1 to avoid conflicts with Julia threading
BLAS.set_num_threads(1)

include("import.jl")
include("kinetic.jl")
include("plots.jl")

using Metaheuristics