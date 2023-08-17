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

"""
In order to simulate the kinetic model, there is information that needs to be provided to the simulator in a specific format:
    1. The initial condition, i.e. the amplitudes of each species/component at the beginning of the simulation. This is given as a vector u0.
    2. Parameter values, i.e. the values of the rate constants. This is given as a vector p.
    3. Time span for which we run the simulation. 

Goals of the GetData function:
    1. Extract all necessary information from the kinetic model defined in the first cell, i.e., the rate constants (k1, k2,...) and components (A, B,...).
    2. Set up vectors u0 and p. 
        a. Create an empty dictionary of the symbolic representations for each component, i.e. key = A, value = []
        b. Do the same for the rate constants, i.e. key = k1, value = []
        b. Map the corresponding values from the input param vector onto the elements in your dictionaries, i.e key = A, value = param[1] = 1
    3. Simulate the kinetic model and convolve the resulting kinetic traces in order to represent the IRF.
    4. Gererate spectral signatures based on these kinetic traces and the real input data.
    5. Generate a simulated 2D data matrix based on the simulated kinetics and spectra.
"""


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


"""
Goals of Objective function:
    1. Read the simulated data matrix from GetData.
    2. Compare simulated data to real data and return the sum of squared differences. 
    3. Return a "residual map", if specified by the user. 
"""

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
function run_optim(obj, bound)
    bboptimize(obj; SearchRange = bound)
end

# generate a residual map
function resMap(param)
    da = GetData(param)
    testData = da[1]
    return heatmap(time, frequencies, (testData - Data), xlabel="Time", ylabel="Wavelength", 
                    title="Residuals Map", colorbar_title="\n \n \n Δ Absorbance", 
                    right_margin=15Plots.mm, left_margin=10Plots.mm, xguidefontsize=10, yguidefontsize=10,
                    c = :thermal)
end
