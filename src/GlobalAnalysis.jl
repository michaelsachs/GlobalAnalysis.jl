module GlobalAnalysis

include("import.jl")
include("kinetic.jl")
include("fits.jl")
include("plots.jl")
include("notebooks.jl")

export DataStruct, SSRMetaData
export importData, importDataVectors, maskData
export getOdeTime, convolveIRF
export generateBounds, generateBoundsBBO, generateBoundsMH
export getSpecies, getParameters, setupVariables
export paramToKin, paramToData, paramToResiduals, paramToResidualsVec
export setupSSRMetaData, projectSSR, paramToSSR, paramToSSRParallel
export getLabels, getParamConfidence, paramToDataCI, printFitResult
export getNumFitVar, calculateAdjR2, calculateChiSq, calculateAICc, calculateBIC
export setup, launch

end #module
