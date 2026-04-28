module GlobalAnalysis

include("import.jl")
include("irf.jl")
include("notebooks.jl")

export importData, getOdeTime, convolveIRF, setup, launch

end #module
