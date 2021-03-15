module GlobalAnalysis

include("TypeDefinitions.jl")
include("DataImport.jl")

#using .TypeDefinitions
using .DataImport

directory = raw"C:\Box Sync\Data\fs-TAS\Matlab data\Fe2O3\400exc\test"

Data = importData(directory)

Data = combineData(Data)

end #module
