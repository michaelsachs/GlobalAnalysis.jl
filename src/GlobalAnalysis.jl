module GlobalAnalysis

export importData

include("TypeDefinitions.jl")
include("DataImport.jl")

directory = raw"C:\Box Sync\Data\fs-TAS\Matlab data\Fe2O3\400exc\test"

Data = importData(directory)


end #module
