using Catalyst
using DifferentialEquations
using NaNStatistics

include("irf.jl")


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
Generates upper and lower limits for fit parameters based on reaction
network `rn` and `limits`, formatted for BlackBoxOptim. 
"""
function generateBoundsBBO(rn, limits)

    lower, upper = generateBounds(rn, limits)

    # convert to array of tuples
    bounds = Array{Tuple{Float64, Float64}}(undef,length(lower))
    for k in eachindex(lower)
        bounds[k] = (lower[k], upper[k])
    end

    return bounds

end

"""
Generates upper and lower limits for fit parameters based on reaction
network `rn` and `limits`, formatted for Metaheuristics. 
"""
function generateBoundsMH(rn, limits)

    lower, upper = generateBounds(rn, limits)
    bounds = [lower'; upper']

    return bounds

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
function simulateData(t, rn, param, limits, Data; ret="res")
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
    sol  = solve(prob, AutoTsit5(Rosenbrock23()))
    kin = transpose(Array(sol))

    # convolve kinetic traces with IRF
    kinConv = convolveIRF(t, kin, μ, σ, tStepParam)

    # normalise kinetics for more consistent output
    # kinConv ./= maximum(kinConv,dims=1)

    # generate spectra based on calculated kinetics
    testSpc = Data / kinConv'
    # assemble data matrix 
    testData = testSpc * kinConv'
    #return testData#, testSpc#, KinMatrix
    if ret == "res"
        return nansum((testData .- Data).^2)
    else
        return testData, testSpc, kinConv
    end

end

"""
Parallel evaluation of `simulateData` for use with Metaheuristics.
`param` is a 2D array, with dimension 1 corresponding to the number
of parallel evaluations, and dimension 2 being the number of fit
parameters.
"""
function simulateDataParallel(t, rn, param, limits, Data)
    fitness = zeros(size(param,1))
    Threads.@threads for n in 1:size(param,1)
        fitness[n] = simulateData(t, rn, param[n,:], limits, Data; ret="res")
    end
    return fitness
end

"""
Generates string array to label species in reaction network `rn`
in plots.
"""
function getLabels(rn)
    permutedims(string.(getSpecies(rn)))
end

