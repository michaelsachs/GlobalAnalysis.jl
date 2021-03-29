
using CSV, Glob
using DataFrames
using StructArrays
using Plots

include("TypeDefinitions.jl")
include("DataImport.jl")

### IMPORT

file = normpath(joinpath(@__DIR__,"..", raw"data\Co3O4_deriv.csv"))
data = CSV.File(file; datarow=2, type=Float64) |> DataFrame!
data = Matrix(data)
wl = data[:,1]
deri1 = data[:,2]
deri2 = data[:,3]

function deri_norm(deri)
    max = maximum(deri)
    min = minimum(deri)
    deri ./= (abs(max) > abs(min) ? abs(max) : abs(min))
end

deri1 = deri_norm(deri1)
deri2 = deri_norm(deri2)

using Interpolations 
itpXYPar = interpolate((XParNaN,YParNaN,), ZParNaN, Gridded(Linear()))


directory = raw"C:\Box Sync\Data\fs-TAS\Julia data\Co3O4\all"
Data = importData(directory)


#### DEFINE PURE SPECTRA

# for now we are not using multiple datasets simultaneously, so we just select one of the datasets
# to be used in the following
data = Data[3]

# y values of the pure z columns. for TAS these are the times of the pure spectra, for SEC the potentials 
# at which pure species would be present 
# this is the simplest way to define pure spectra. you'll probably want to add the option to import 
# spectra from a csv file as you don't necessarily have times/potentials at which only one species is
# present. These then most likely will need to be interpolated to match the spectral sampling of the SEC dataset
pureY = [0.4, 0.9] #[0.7,6000]

pureZ = Array{Union{Float64,Missing}}(undef,length(data.x),2)
for k in 1:length(pureY)
    pureZIdx = argmin(abs.(data.y .- pureY[k]))
    pureZ[:,k] = data.z[:,pureZIdx]
end



#### OPTIMISATION 
#=
Below two different ways are given to calculate kinetics of using the pure spectra. Both independently
fit spectra at each time/potential using a linear combination of the pure component spectra. The 
NLopt way does this sequentially for each spectrum, whereas BlackBoxOptim fits all spectra 
simultaneously but still in an independent manner. I would recomment NLopt as it seems much faster 
=#



#### OPTIMISATION 1: NLOPT

using NLopt
using LinearAlgebra: Diagonal

# this function calculates residuals for the optimisation. Loss functions in NLopt have a gradient as 
# mandatory input; we are not using one here, but the function still needs to be able to 
# accept grad as a parameter
function lossFun(param, grad)
    testZ = pureZ * Diagonal(param)
    #coalesce sets missings to zero for fitting purposes
    #sum(abs2,x-y) calculates sum of squares
    loss = sum(abs2, coalesce.(sum(testZ,dims=2),0) .- coalesce.(currentZ,0))
end

# set up boundaries and start parameters for an intial guess
lower = fill(0.0, length(pureY))
upper = fill(10.0, length(pureY))
startParam = fill(0.5, length(pureY))

# pre-allocate array to collect fit result
paramOpt = Array{Float64}(undef,length(data.y),length(pureY))
# time the loop execution just to keep track of performance
@time begin
for ny in 1:length(data.y)
    # y slice in z matrix (e.g. current spectrum)
    currentZ = data.z[:,ny]

    # set up NLopt parameters
    # COBYLA is a gradient free optimisation algorithm, hence grad is unused above
    opt = Opt(:LN_COBYLA, length(startParam)) 
    lower_bounds!(opt,lower)
    upper_bounds!(opt,upper)
    min_objective!(opt, lossFun)
    xtol_rel!(opt, 1e-8)
    maxeval!(opt, 1000)

    # run optimisation 
    (NLOminf,NLOminx,NLOret) = NLopt.optimize(opt, startParam)
    @show NLOret
    # write fitted paramers into paramOpt
    paramOpt[ny,:] = NLOminx

end
end #time end


# plot the result
# boolean array to remove all data ≤ 0, otherwise error thrown when making x axis logarithmic
# the dot . is a broadcast operator, meaning that I want to apply > 0 to all elements in data.y
yBool = data.y .> 0
# make a plot like this
p1 = plot(data.x,pureZ, linewidth=2, xlabel="Wavelength(nm)", ylabel = "ΔA", title="Pure spectra")
# and another one
p2 = plot(data.y[yBool], paramOpt[yBool,:],xscale=:log10,linewidth=2,xlabel="Time (ps)", 
    ylabel = "Amplitude (a.u.)", title="Kinetics from spectral model")

# assemble 2 panel plot like this
plot(p1,p2,layout=(2,1))
# modify a plot like this (! in julia functions indicates that they modify the input in-place)
plot!(legend=false, titlefontsize=12, size=(600,600))






#### OPTIMISATION 2: BBO

using BlackBoxOptim

# BBO takes 1D array of tuples as input. However we need a 2D array if we 
# have 2 or more pure spectra (i.e. always), hence the parameter vector is formulated 
# as a 1D array here and then reshaped in the loss function 

# bounds here need to be set up as tuples
bounds = fill((0.0,1.0),length(data.y)*length(pureY))

function lossFun(param)
    param = reshape(param, (length(data.y),length(pureY))) 
    #matrix multiplication to create guess matrix to compare to actual data
    testZ = pureZ * transpose(param)
    #set missings to zero for fitting purposes
    loss = sum(abs2, coalesce.(testZ,0) .- coalesce.(data.z,0))
end

# upper time limit in seconds for the fit. Probably needs several minutes for accurate results
MaxTimeBBO = 60
resultBBO = bboptimize(lossFun; SearchRange = bounds, NThreads=Threads.nthreads()-1, 
    MaxTime=MaxTimeBBO, Method = :adaptive_de_rand_1_bin_radiuslimited, maxiters=50)

# extract optmimised paraetmers 
paramOpt = best_candidate(resultBBO)
paramOpt = reshape(paramOpt, (length(data.y),length(pureY))) 

yBool = data.y .> 0
plot(data.y[yBool], paramOpt[yBool,:],xscale=:log10)




