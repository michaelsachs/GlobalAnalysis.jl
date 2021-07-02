include("TypeDefinitions.jl")
include("DataImport.jl")
include("IRFConvolution.jl")


directory = raw"C:\Box Sync\Data\fs-TAS\Julia data\Co3O4\all"

Data = importData(directory; miss="NaN")
Data = importData(directory)
nDat = 3

X = coalesce.(Data[nDat].x,0)
Y = coalesce.(Data[nDat].y,0)
#convert to ΔT/T
Z = coalesce.(10 .^ -Data[nDat].z .- 1,0)

function data2txt(data, io)
    for k in 1:size(data,2)
        for n in 1:size(data,1)
            print(io, "$(data[n,k]) ")
        end
        print(io, "\n")
    end

end

filename = "$(Data[nDat].name)_$(Data[nDat].val[1])$(Data[nDat].unit[1])"

try
    rm("$(filename).ana")
catch
end

open("$(filename).ana", "w") do io

    println(io, "%FILENAME=$(filename)")
    println(io, "%DATATYPE=TAVIS")
    println(io, "%TIMESCALE=ps")

    print(io, "%TIMELIST=")
    data2txt(Y, io)

    print(io, "%WAVELENGTHLIST=")
    data2txt(X, io)

    print(io, "%INTENSITYMATRIX=\n")
    data2txt(Z, io)

end