#imports data from .csv, automatically combines datasets based on used parameters and x/y ranges 
#over which they contain data (e.g. different or overlapping probe ranges). To this end data is 
#interpolated if required. Result is datastruct which contains one element for each unique set 
#of parameters (e.g. fluence, temperature)

include("TypeDefinitions.jl")

using CSV, Glob
using DataFrames
using StructArrays
using Statistics
using Interpolations


function dataOverlap(XPar)
    #bool for non-overlapping parts of vector
    overlapBoolAr = []

    for nPar in 1:length(XPar)
        #create view of other datasets than current one
        XParBool = trues(length(XPar))
        XParBool[nPar] = false
        XParOther = @view XPar[XParBool]
        
        #bool for non-overlapping parts of vector
        overlapBool = trues(length(XPar[nPar]))
        for nParOther in 1:length(XParOther)
            XParOtherBool = (XPar[nPar] .< XParOther[nParOther][1]) .| (XPar[nPar] .> XParOther[nParOther][end])
            overlapBool = overlapBool .& XParOtherBool
        end

        push!(overlapBoolAr, .!overlapBool)

    end
    return overlapBoolAr
end

# determine vector to interpolate other datasets onto 
# returns vector with highest number of occurences along dim
# take into account spacing??
function determineItpVector(XPar, XParItpBoolAr)
    XParUnique = []
    XParUniqueNum = []
    for nn in 1:length(XPar)

        XParItp = XPar[nn][XParItpBoolAr[nn]]
        # check if vector is already part of XParUnique
        XParItpBool = XParUnique .== [XPar[nn]]
        if XParItp == []
            continue
        elseif any(x->x==true, XParItpBool) == true
            XParUniqueNum[XParItpBool] .+= length(XParItp)
        else
            push!(XParUnique, XParItp)
            push!(XParUniqueNum, length(XParItp))
        end
    end
    XParCommon = XParUnique[sortperm(XParUniqueNum)][end]
end

# returns bool of length A to identify position of a in A
function findVec(A,a)
    for n in 0:length(A)-length(a)
        bool = [falses(n); trues(length(a)); falses(length(A)-length(a)-n)]
        if A[bool] == a
            return bool
        end
    end
end

# imports all .csv files in chosen directory
# first column in csv is x data, first coumn is y data, [1,1] is not used, rest is z data
function importData(directory)
    files = glob("*.csv", directory)
    fileNames = readdir(directory)
    # remove elements which do not end in .csv
    filter!(x->x[end-3:end]==".csv", fileNames)

    X = []
    Y = []
    Z = []
    param = []
    for nFile in 1:length(files)

        #import; all entries which cannot be converted to Float64 are imported as missing
        DF = CSV.File(files[nFile]; datarow=1, type=Float64) |> DataFrame!
        #replace NaN with missing
        allowmissing!(DF)
        [replace!(col, NaN=>missing) for col = eachcol(DF)]
        #convert DataFrame to array. At this point it should only contain Float64 and missing
        DF = Matrix(DF)

        #count missings in each col
        NumMissingsCols = [sum(ismissing.(col)) for col = eachcol(DF)]
        #count missings in each row
        NumMissingsRows = [sum(ismissing.(row)) for row = eachrow(DF)]
        #remove cols/rows with >90% missings
        DF = DF[NumMissingsRows .< (size(DF,2)*0.90), NumMissingsCols .< (size(DF,1)*0.90)]

        #separate imported data into x, y, z 
        x = DF[2:end, 1]
        y = DF[1, 2:end]
        z = DF[2:end, 2:end]
        
        
        # average x/y entries with same value 
        xUnique = unique(x)
        if x ≠ xUnique
            zxUnique = Array{Union{Float64,Missing}}(undef,length(xUnique),length(y))
            for nx in 1:length(xUnique)
                xBool = x .== xUnique[nx]
                zxUnique[nx,:] = mean(z[xBool,:], dims=1)
            end
        else
            zxUnique = @view z[:,:]
        end

        yUnique = unique(y)
        if y ≠ yUnique
            zUnique = Array{Union{Float64,Missing}}(undef,length(xUnique),length(yUnique))
            for ny in 1:length(yUnique)
                yBool = y .== yUnique[ny]
                zUnique[:,ny] = mean(zxUnique[:,yBool], dims=2)
            end
        else
            zUnique = @view z[:,:]
        end
       
        
        # extract variable data (e.g. fluence; temperature)
        # divide by separators "__"
        # add m in case string starts with "__"
        strArray = split("m"*fileNames[nFile], "__")
        # no separator
        if length(strArray) == 1
            name = fileNames[nFile]
            vars = []
            vals = []
            units = []
        # one separator
        elseif length(strArray) == 2
            error("\"$(fileNames[nFile])\" contains one \"__\" separator. Must contain ≥2 or none.")
        # two separators
        elseif length(strArray) == 3
            name = strArray[2]
            vars = []
            vals = []
            units = []
        # at least 2 separators
        elseif length(strArray) ≥ 4
            strArray = strArray[2:end-1]
            name = strArray[1]
            strArrayAllVar = strArray[2:end]
            vars = Array{String}(undef, length(strArrayAllVar))
            vals = Array{Float64}(undef, length(strArrayAllVar))
            units = Array{String}(undef, length(strArrayAllVar))
            for numVar in 1:length(strArrayAllVar)
                #divide into var and unit
                strArraySingleVar = split(strArrayAllVar[numVar], "_")
                vars[numVar] = strArraySingleVar[1]
                try
                    vals[numVar] = parse(Float64, strArraySingleVar[2])
                catch err
                    error("\"$(strArraySingleVar[2])\" in \"$(fileNames[nFile])\" must not contain non-numeric characters")
                end

                if length(strArraySingleVar) == 2
                    units[numVar] = ""
                elseif length(strArraySingleVar) == 3
                    units[numVar] = strArraySingleVar[3]
                else
                error("\"$(strArray[numVar])\" must not contain more than two \"_\" separators")
                end
            end
            #sort alphabetically by var
            ids = sortperm(vars)
            vars = vars[ids]
            vals = vals[ids]
            units = units[ids]
        end

        push!(X, xUnique)
        push!(Y, yUnique)
        push!(Z, zUnique)
        push!(param, [name,vars,vals,units])

    end
    
    paramUnique = unique(param)

    ### combine datasets which have same parameters
    Data = StructArray{DataStruct}(undef,length(paramUnique))
    for nAllParam in 1:length(paramUnique)
        paramBool = [paramUnique[nAllParam]] .== param

        XPar = X[paramBool]
        YPar = Y[paramBool]
        ZPar = Z[paramBool]


        #delete identical datasets

        #need to check overlaps first
        #divide into regions of 0/1/2D overlap

        #divide data for same parameter into data blocks
        #1) no overlap: append blocks
        #2) 1D overlap: interpolate overlapping blocks
            #use vec with which occurs most times?
        #3) 2D overlap??
            #average if points are exacly same 
            #otherwise pick one vector for interpolation 

        # interpolate if XPar or YPar elements are not all same
        if (all(x->x==XPar[1], XPar) == false) | (all(y->y==YPar[1], YPar) == false)
            #bool arrays where true indicates that interpolation is needed
            XParOverlapBoolAr = dataOverlap(XPar)
            YParOverlapBoolAr = dataOverlap(YPar)

            #vectors to interpolate onto ####combine both functions 
            XParCommon = determineItpVector(XPar, XParOverlapBoolAr)
            YParCommon = determineItpVector(YPar, YParOverlapBoolAr)
            #####temp
            #XParCommon = [XParCommon; collect(1000:10:1100)]
            #YParCommon = [YParCommon; collect(6000:1000:7000)]


            
            # allows for common itp vectors to contain arbitrary elements without being bound to array dimensions (future)
            #could be condensed into function....
            #need to add equal sign to one of the 2
            XParItpBoolAr = []
            YParItpBoolAr = []
            XParCommonBoolAr = []
            YParCommonBoolAr = []
            for bb in 1:length(XPar)
                #interpolate if overlap with other datasets OR within common itp dimensions, BUT only within dimensions of
                #current param dataset
                XParItpBool = (XPar[bb] .> XParCommon[1]) .& (XPar[bb] .< XParCommon[end]) .| XParOverlapBoolAr[bb] ###emprty
                YParItpBool = (YPar[bb] .> YParCommon[1]) .& (YPar[bb] .< YParCommon[end]) .| YParOverlapBoolAr[bb]
                push!(XParItpBoolAr, XParItpBool)
                push!(YParItpBoolAr, YParItpBool)

                XParCommonBool = (XParCommon .> XPar[bb][1]) .& (XParCommon .< XPar[bb][end]) ###empty
                YParCommonBool = (YParCommon .> YPar[bb][1]) .& (YParCommon .< YPar[bb][end]) 
                push!(XParCommonBoolAr, XParCommonBool)
                push!(YParCommonBoolAr, YParCommonBool)
            end

            #if empty: make overlap all full and commonBool all empty

            
            # interpolation onto YParCommon. Note that missing is not currently supported by Interpolations, therefore missings are 
            # for now converted to NaN and then back to missing after interpolation 
            ZParXInterpolated = []
            ZParYInterpolated = []
            ZParXYInterpolated = []
            for kk in 1:length(XPar) #####no overlap case?? Xparnan empty. check zparall (should be appended)
                XParNaN = coalesce.(XPar[kk][XParOverlapBoolAr[kk]], NaN) #chnage way these emtpy ones are stored?
                YParNaN = coalesce.(YPar[kk][YParOverlapBoolAr[kk]], NaN)
                ZParNaN = coalesce.(ZPar[kk][XParOverlapBoolAr[kk],YParOverlapBoolAr[kk]], NaN)
                XParCommonNaN = coalesce.(XParCommon[XParCommonBoolAr[kk]],NaN) 
                #emty means there are no common ones; all should be appended later

                YParCommonNaN = coalesce.(YParCommon[YParCommonBoolAr[kk]],NaN)
                #itpYPar = interpolate((XPar[kk],YPar[kk],), ZPar[kk], Gridded(Linear()))
                #replace missing by NaN as long as it is not supported by Interpolations
                itpXYPar = interpolate((XParNaN,YParNaN,), ZParNaN, Gridded(Linear()))
                #extrapolate to NaN outside array dimensions. Replace by missing once supported
                itpXYPar = extrapolate(itpXYPar, NaN)
                # separate interpolated vectors for X and Y to avoid returning only their overlap
                ZParXinterp = itpXYPar(XParCommonNaN, YParNaN)
                ZParXinterp = replace(ZParXinterp, NaN=>missing)
                push!(ZParXInterpolated, ZParXinterp)

                ZParYinterp = itpXYPar(XParNaN, YParCommonNaN)
                ZParYinterp = replace(ZParYinterp, NaN=>missing)
                push!(ZParYInterpolated, ZParYinterp)

                ZParXYinterp = itpXYPar(XParCommonNaN, YParCommonNaN)
                ZParXYinterp = replace(ZParXYinterp, NaN=>missing)
                push!(ZParXYInterpolated, ZParXYinterp)
            end

            itpXYPar = interpolate((XParNaN,YParNaN,), ZParNaN, Gridded(Linear()))


            # interpolation onto YParCommon. Note that missing is not currently supported by Interpolations, therefore missings are 
            # for now converted to NaN and then back to missing after interpolation 
            ZParXInterpolated = []
            ZParYInterpolated = []
            ZParXYInterpolated = []
            for kk in 1:length(XPar)
                XParNaN = coalesce.(XPar[kk], NaN)
                YParNaN = coalesce.(YPar[kk], NaN)
                ZParNaN = coalesce.(ZPar[kk], NaN)
                XParCommonNaN = coalesce.(XParCommon[XParCommonBoolAr[kk]],NaN)
                YParCommonNaN = coalesce.(YParCommon[YParCommonBoolAr[kk]],NaN)
                #itpYPar = interpolate((XPar[kk],YPar[kk],), ZPar[kk], Gridded(Linear()))
                #replace missing by NaN as long as it is not supported by Interpolations
                itpXYPar = interpolate((XParNaN,YParNaN,), ZParNaN, Gridded(Linear()))
                #extrapolate to NaN outside array dimensions. Replace by missing once supported
                itpXYPar = extrapolate(itpXYPar, NaN)
                # separate interpolated vectors for X and Y to avoid returning only their overlap
                ZParXinterp = itpXYPar(XParCommonNaN, YParNaN)
                ZParXinterp = replace(ZParXinterp, NaN=>missing)
                push!(ZParXInterpolated, ZParXinterp)

                ZParYinterp = itpXYPar(XParNaN, YParCommonNaN)
                ZParYinterp = replace(ZParYinterp, NaN=>missing)
                push!(ZParYInterpolated, ZParYinterp)

                ZParXYinterp = itpXYPar(XParCommonNaN, YParCommonNaN)
                ZParXYinterp = replace(ZParXYinterp, NaN=>missing)
                push!(ZParXYInterpolated, ZParXYinterp)
            end

            
            # assemble new ZPar arrays from interpolated and non-interpolated parts of each dataset
            # and sort arrays 
            XParNew = []
            YParNew = []
            ZParNew = []
            for kk in 1:length(XPar)
                XParUnsorted = [XParCommon[XParCommonBoolAr[kk]]; XPar[kk][.!XParOverlapBoolAr[kk]]]
                YParUnsorted = [YParCommon[YParCommonBoolAr[kk]]; YPar[kk][.!YParOverlapBoolAr[kk]]]
                XParSortBool = sortperm(XParUnsorted)
                YParSortBool = sortperm(YParUnsorted)
                push!(XParNew, XParUnsorted[XParSortBool])
                push!(YParNew, YParUnsorted[YParSortBool])
                
                XYItp_YItp = [ZParXYInterpolated[kk]; ZParYInterpolated[kk][.!XParOverlapBoolAr[kk],:]]
                XItp_NoItp = [ZParXInterpolated[kk][:,.!YParOverlapBoolAr[kk]]; ZPar[kk][.!XParOverlapBoolAr[kk],.!YParOverlapBoolAr[kk]]]
                ItpFinal = hcat(XYItp_YItp,XItp_NoItp)
                ItpFinalSorted = ItpFinal[XParSortBool, YParSortBool]
                push!(ZParNew, ItpFinalSorted)
            end

            ### combine parameter datasets into one 
            XParAll = sort!(unique!(vcat(XParNew...)))
            YParAll = sort!(unique!(vcat(YParNew...)))
            ZParAll = zeros(Union{Float64,Missing}, length(XParAll), length(YParAll))
            #array to keep track of how many datasets contribute to each datapoint
            numZParAll = zeros(Union{Int64,Missing}, length(XParAll), length(YParAll))

            for öö in 1:length(XPar)
                XParBool = findVec(XParAll,XParNew[öö])
                YParBool = findVec(YParAll,YParNew[öö])

                ZParAll[XParBool,YParBool] .+= ZParNew[öö]
                numZParAll[XParBool,YParBool] .+= 1
            end
            # in the next step this turns all ZParAll points with no data into missings
            replace!(numZParAll, 0=>missing)
            # new array where all points with contributions from multiple datasets are
            # averaged
            ZParAll ./= numZParAll

        # all X and Y are the same, hence no interpolation required
        else
            XParAll = XPar[1]
            YParAll = YPar[1]

            ZParAll = zeros(Union{Float64,Missing}, length(XParAll), length(YParAll))
            for pp in 1:length(XPar)
                ZParAll .+= ZPar[pp]
            end
            ZParAll ./= length(XPar)
        end
        
        #store combined data in new DataStruct
        #if length(paramUnique[nAllParam][3]) == 1
        #Data[nAllParam] = deepcopy(DataStruct(paramUnique[nAllParam][1], XParAll, YParAll, ZParAll, 
        #    paramUnique[nAllParam][2], paramUnique[nAllParam][3], paramUnique[nAllParam][4], 
        #    [fileNames[paramBool]]))
        #else
        Data[nAllParam] = deepcopy(DataStruct(String(paramUnique[nAllParam][1]), XParAll, YParAll, ZParAll, 
            paramUnique[nAllParam][2], paramUnique[nAllParam][3], paramUnique[nAllParam][4], 
            fileNames[paramBool]))
        #end

    end
    return Data
end

#make X and Y non missings?

