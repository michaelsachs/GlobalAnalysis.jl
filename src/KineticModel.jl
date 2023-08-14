using Catalyst, Plots, DifferentialEquations
using BlackBoxOptim
using CSV
using DataFrames
using DataStructures
using Statistics
using NaNStatistics

include("IRFConvolution.jl")
include("DataImport.jl")
include("TypeDefinitions.jl")
####### ONLY EDIT BETWEEN THESE LINES #######

RealData = importData("/Users/jessicaflowers/Desktop/JULIA/data/CdO/"; miss="NaN")
RealData = RealData[1]

# define a kinetic model
rs = @reaction_network begin
    k1, A --> B
    k2, B --> 0
end

# define bounds for the parameter optimization 
state_lower = [1, 0] # A, B, C
state_upper = [1, 0] 

IRF_lower = [0.1, 0.01]
IRF_upper = [0.2, 0.08]

rate_const_lower = [1, 0.0001]
rate_const_upper = [3, 0.0003]


####### ONLY EDIT BETWEEN THESE LINES #######

# some pre-processing of the input csv file
frequencies = RealData.x
# frequencies = vcat(frequencies...)
time = RealData.y
# time = vcat(time...)
Data = RealData.z
# Data = (vcat(Data...))
# Data = coalesce.(Data, NaN)


# create bounds array
lower = vcat(state_lower, IRF_lower, rate_const_lower)
upper = vcat(state_upper, IRF_upper, rate_const_upper)

bounds = Array{Tuple{Float64, Float64}}(undef,length(lower))
for k in 1:length(lower)
    bounds[k] = (lower[k], upper[k])
end

# function to produce data with "guessed" parameters
function GetData(param) # A, B, ..., μ, σ, k1, k2, k3
    # Extract the number of states/components and rate constants from the kinetic model
    num_states = length(states(rs)) # number of components
    num_k = length(Catalyst.parameters(rs)) # number of rate constants
    States = states(rs) # get components from kinetic model, A(t), B(t), C(t) ...

    u0_ = OrderedDict{Symbol, Float64}() # initialize the component amplitude vector
    for (i, state) in enumerate(States) # loop through the components A, B ... extracted from kinetic model
        name = string(state)
        name = match(r"(\w+)", name).captures[1]  # Extract alphanumeric part of the name, i.e. A rather than A(t)
        u0_[Symbol(name)] = param[i] # assign the corresponding value from optimized param vector
    end
    u0 = Pair{Symbol, Float64}[]
    for (key, value) in pairs(u0_)
        push!(u0, Pair(key, Float64(value))) 
    end
    u0 = collect(u0) # collect initial value vector for mapping onto parameters 

    num_p = num_states + num_k + 2 # the number of parameters to optimize is the number of states + rate constants + μ + σ
    
    k_start = num_states + 3 # rate constants in parameter vector start after the states and μ and σ, i.e. num_states + 3
    unsorted_k = ((param[k_start:num_p])) # rate constants, unsorted

    # get the rate constants from the kinetic model, k1, k2, k3
    Ks = Catalyst.parameters(rs)
    p_ = OrderedDict{Symbol, Float64}()

    for (i, K) in enumerate(Ks) # loop through rate constants extracted from kinetic model
        rate_const = string(K)
        p_[Symbol(rate_const)] = unsorted_k[i]
    end

    p = Pair{Symbol, Float64}[]
    for (key, value) in pairs(p_)
        push!(p, Pair(key, Float64(value)))
    end
    collect(p) # collect rate constants into vector

    # IRF parameters
    # first IRF parameter is after the component concentrations, i.e. after num_states
    μ=param[num_states+1]
    σ=param[num_states+2]

    TimeODE, TimeIRFPos, IRFWindowCut, TimeODEBool = defineIRFTime(time, μ, σ)
    tspan = [minimum(TimeODE),maximum(TimeODE)] 

    prob = ODEProblem(rs, u0, tspan, p; saveat=TimeODE)
    sol  = solve(prob, Tsit5())
    sol_array = Array(sol)

    # Perform convolution for each component
    ConvKin = Vector{Array{Float64,1}}(undef, num_states)
    for i in 1:num_states # loop through components
        component = sol_array[i,:]
        ConvKin[i] = convolveIRF(time, component, μ, σ, TimeIRFPos, IRFWindowCut, TimeODEBool)
    end

    KinMatrix = reduce(vcat, ConvKin')  # kinetic matrix
    TestSpec = Data / KinMatrix  # spectral signatures matrix
    testData = TestSpec * KinMatrix # test matrix from kinetics and signatures
    return testData, TestSpec, KinMatrix
end

# objective function which returns the sum of squared differences bewteen 
# the produced data with guessed parameters (GetData) and the real data matrix
function Objective(param; output="res")
    da = GetData(param)
    testData = da[1]
    if output == "res"
        return nansum((testData .- Data).^2) #sum(abs2, testData .- Data)
    elseif output == "map"
        return heatmap(time, frequencies, (testData - Data), xlabel="Time", ylabel="Wavelength", 
                    title="Residuals Map", colorbar_title="\n \n \n Δ Absorbance", 
                    right_margin=15Plots.mm, left_margin=10Plots.mm, xguidefontsize=10, yguidefontsize=10,
                    c = :thermal)
    end
end

# minimize objective function to optimize parameters
# OP = bboptimize(Objective; SearchRange = bounds)
RecoveredParam = best_candidate(OP)

RecoveredData = GetData(RecoveredParam)
DataMatrix = RecoveredData[1] # data matrix
RecoveredSpec = RecoveredData[2] # spectral signatures
RecoveredKinetics = RecoveredData[3] # kinetic traces

# kinetics
tpos = time.>0
plot(time[tpos], RecoveredKinetics[1,:][tpos], xscale=:log10,
     title="Recovered Kinetics", xlabel="log(Time)", legend=:topleft, linewidth = 2)
plot!(time[tpos], RecoveredKinetics[2,:][tpos], xscale=:log10, linewidth = 2)
# plot!(time[tpos], RecoveredKinetics[3,:][tpos], xscale=:log10, linewidth = 2)



# spectral signatures
plot(frequencies, RecoveredSpec, xlabel="Wavelength", linewidth = 2) 

plot(frequencies, RecoveredSpec, title="Recovered Spectral Signatures", xlabel="Wavelength", linewidth = 2) 

# 3D data 
surface(time, frequencies, DataMatrix, title="Recovered Data", xlabel="Time", ylabel="         Wavelength", xguidefontsize=10, yguidefontsize=10, xrotation = 45, right_margin=15Plots.mm, colorbar_title="\n\n\n\nΔ Absorbance")
surface(time, frequencies, Data, title="Real Data", xlabel="Time", ylabel="         Wavelength", xguidefontsize=10, yguidefontsize=10, xrotation = 45, right_margin=15Plots.mm, colorbar_title="\n\n\n\nΔ Absorbance")

# residual maps 
map = Objective(RecoveredParam; output="map")

real_heatmap = heatmap(time, frequencies, Data, 
    xlabel="Time", ylabel="Wavelength", title="Real Data", 
    colorbar_title="\n\nΔ Absorbance", right_margin=15Plots.mm, left_margin=10Plots.mm, 
    xguidefontsize=10, yguidefontsize=10,
    c = :thermal)
optim_heatmap = heatmap(time, frequencies, DataMatrix, 
    xlabel="Time", ylabel="Wavelength", title="Recovered Data" ,
    colorbar_title="\n\nΔ Absorbance", right_margin=15Plots.mm, left_margin=10Plots.mm, 
    xguidefontsize=10, yguidefontsize=10,
    c = :thermal)

residual_val = Objective(RecoveredParam; output="res")

# # Normalize the spectral signatures and kinetics

# kn1 = RecoveredKinetics[1,:]
# kn2 = RecoveredKinetics[2,:]
# kn3 = RecoveredKinetics[3,:]

# max1 = maximum(kn1)
# max2 = maximum(kn2)
# max3 = maximum(kn3)


# kn1 = kn1 / max1
# kn2 = kn2 / max2
# kn3 = kn3 / max3

# plot(kn1)
# plot!(kn2)
# plot!(kn3)

# spc1 = RecoveredSpec[:,1]
# spc2 = RecoveredSpec[:,2]
# spc3 = RecoveredSpec[:,3]

# Max1 = maximum(filter(!isnan, spc1))
# Max2 = maximum(filter(!isnan, spc2))
# Max3 = maximum(filter(!isnan, spc3))

# spc1 = spc1 / Max1
# spc2 = spc2 / Max2
# spc3 = spc3 / Max3

# plot(spc1)
# plot!(spc2)
# plot!(spc3)
