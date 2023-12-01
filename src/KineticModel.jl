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

"""
Generates upper and lower limits for fit parameters based on reaction
network `rn` and `limits`. 
"""
function generateBounds(rn, limits)
    syms = [getSpecies(rn); getParameters(rn); :μ; :σ]

    count = 0
    lower = Float64[]
    upper = Float64[]
    for sym in syms
        if sym ∈ keys(limits)
            val = limits[sym]
            # fit parameter is indicated by range
            if length(val) == 2
                count += 1
                push!(lower,minimum(val))
                push!(upper,maximum(val))
            end
        end
    end
    return lower, upper
end


"""
Assembles parameter dictionaries for ODE solver, drawing values from (i) fit 
parameter vector `fitParam` for fit variables, (ii) the `limits` dictionary if 
the variable has a fixed value, or (iii) a `default` value for variables which 
are neither fitted nor supplied in `limits`. `count` keeps track of how many 
fit parameters have been assigned.
"""
function gatherParams(syms, fitParam, limits, default, count)
    p = Dict{Symbol,Float64}()
    for sym in syms
        if sym ∈ keys(limits)
            val = limits[sym]
            # fixed parameter
            if length(val) == 1
                if val isa Number
                    p[sym] = val
                elseif val isa Vector
                    p[sym] = val[1]
                end
            # fit parameter
            elseif length(val) == 2
                count += 1
                p[sym] = fitParam[count]
            else
                error("Bounds value for $(sym) must be
                (a) a number if $(sym) is a fixed parameter, or
                (b) a vector containing two numbers if $(sym) is a fit parameter")
            end
        # default parameter
        else
            p[sym] = default
        end
    end
    return p, count
end

"""
Returns species as a vector of symbols.
"""
function getSpecies(rn)
    species = Catalyst.states(rn)
    # species as strings
    speciesStr = map(x->match(r"(\w+)",x).captures[1], string.(species))
    # species as symbols
    return Symbol.(speciesStr)
end

"""
Returns parameters as a vector of symbols.
"""
function getParameters(rn)
    rateConst = Catalyst.parameters(rn)
    # rate constants as symbols
    return Symbol.(rateConst)
end


"""
Generates test data by calculating kinetic traces based on time vector
`t`, reaction network `rn`, and parameter vector `param`. `param` contains 
parameters in order (1) initial state populations, (2) rate constants, 
(3) IRF parameters.
"""
function simulateData(t, rn, param, limits, Data)
    species = getSpecies(rn)
    # get populations at t = 0 
    u0, count = gatherParams(species, param, limits, 1, 0)

    rateConst = getParameters(rn)
    # get rate constants
    ks, count = gatherParams(rateConst, param, limits, NaN, count)

    # get Gaussian IRF parameters
    irfParam, count = gatherParams([:μ,:σ], param, limits, 0, count)
    μ = irfParam[:μ]
    σ = irfParam[:σ]

    # make sure that all fit parameters have been distributed
    @assert length(param) == count

    # assemble time vector for ODE solver
    tStepParam = getOdeTime(t, μ, σ)
    tOde = tStepParam[1]
    # time span for ODE solver
    tspan = [minimum(tOde), maximum(tOde)] 

    # set up and solve ODEs
    prob = ODEProblem(rn, u0, tspan, ks; saveat=tOde)
    sol  = solve(prob, ROS3P())
    kin = transpose(Array(sol))

    # convolve kinetic traces with IRF
    kinConv = convolveIRF(t, kin, μ, σ, tStepParam)

    # generate spectra based on calculated kinetics
    testSpc = Data / kinConv'
    # assemble data matrix 
    testData = testSpc * kinConv'
    return testData#, testSpc#, KinMatrix
end


"""
Goals of Objective function:
    1. Read the simulated data matrix from GetData.
    2. Compare simulated data to real data and return the sum of squared differences. 
    3. Return a "residual map", if specified by the user. 
"""

function Objective(param; output="res")
    da = simulateData(param)
    testData = da#[1]
    if output == "res"
        return nansum((testData .- Data).^2) #sum(abs2, testData .- Data)
    elseif output == "map"
        return heatmap(time, wavelength, (testData - Data), xlabel="Time", ylabel="Wavelength", 
                    title="Residuals Map", colorbar_title="\n \n \n Δ Absorbance", 
                    right_margin=15Plots.mm, left_margin=10Plots.mm, xguidefontsize=10, yguidefontsize=10,
                    c = :thermal)
    end
end

# minimize objective function to optimize parameters
function run_optim(obj, bound)
    bboptimize(obj; SearchRange = bound)
end
