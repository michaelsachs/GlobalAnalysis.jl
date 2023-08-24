module GlobalAnalysis

#export importData, defineIRFTime, convolveIRF, maskData, jupyterNotebook, jupyterLab
export jupyterNotebook, jupyterLab

#include("TypeDefinitions.jl")
#include("DataImport.jl")
#include("IRFConvolution.jl")
include("notebooks.jl")


end #module
