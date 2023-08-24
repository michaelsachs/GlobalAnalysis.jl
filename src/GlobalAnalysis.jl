module GlobalAnalysis

export importData, defineIRFTime, convolveIRF, maskData, jupyterNotebook, jupyterLab

include("TypeDefinitions.jl")
include("DataImport.jl")
include("IRFConvolution.jl")
include("notebooks.jl")

#directory = raw"C:\Box Sync\Data\fs-TAS\Matlab data\Fe2O3\400exc\test"

#Data = importData(directory; miss="NaN")
#Data = importData(directory)


end #module
