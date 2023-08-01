using Catalyst, Plots, DifferentialEquations
using BlackBoxOptim
using CSV
using DataFrames
# using DataStructures
include("IRFConvolution.jl")
include("DataImport.jl")
include("TypeDefinitions.jl")
####### ONLY EDIT BETWEEN THESE LINES #######

# da = CSV.read("/Users/jessicaflowers/Desktop/JULIA/data/Fe2O3_0.20.csv", DataFrame)
testData = importData("/Users/jessicaflowers/Desktop/JULIA/data/")
# define a kinetic model
rs = @reaction_network begin
    k1, A --> 0
    k2, 2B --> 0
    k3, C --> 0
end

# define bounds for the parameter optimization 
state_lower = [1, 1, 1] # A, B, C
state_upper = [1, 1, 1] # for first order, set the first value to 1 and the rest to 0

IRF_lower = [0.3, 0.04]
IRF_upper = [0.5, 0.06]

rate_const_lower = [0.01, 0.001, 0.00001]
rate_const_upper = [0.9, 0.05, 0.00005]

####### ONLY EDIT BETWEEN THESE LINES #######

# some pre-processing of the input csv file
# convert dataframe into matrix
dat = Matrix(da)
frequencies = dat[:,1]
data = dat[:, 1:end .≠ 1] # get rid of 1st column as this is just the freq range, not real data


function parse_time_column(column_name::String)
    parsed_string = split(column_name, "_")[1]  # Take only the part before "_"
    return parse(Float64, parsed_string)
end

# Extract the time stamps and convert them into a vector of Float64
col_names = names(da)
time_stamp = [parse_time_column(col_name) for col_name in col_names[2:end]]

# create bounds array
lower = vcat(state_lower, IRF_lower, rate_const_lower)
upper = vcat(state_upper, IRF_upper, rate_const_upper)

bounds = Array{Tuple{Float64, Float64}}(undef,length(lower))
for k in 1:length(lower)
    bounds[k] = (lower[k], upper[k])
end

TimeDiff = diff(time_stamp)
IRFStep = minimum(filter(t -> t > 0, TimeDiff))
IRFWindow = 30*σ
TimeIRFPos = collect(IRFStep/2:IRFStep:IRFWindow)
IRFWindowCut = IRFWindow*0.5
TimeODE = sort(unique([TimeIRFPos; time_stamp[time_stamp .> IRFWindowCut] .- μ]))
TimeODEBool = [in.(TimeODE, Ref(Set(TimeIRFPos))), in.(TimeODE, Ref(Set(time_stamp .- μ)))]

tspan = [minimum(TimeODE),maximum(TimeODE)] 

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

    prob = ODEProblem(rs, u0, tspan, p; saveat=TimeODE)
    sol  = solve(prob, Tsit5())
    sol_array = Array(sol)
    # kinODE = vcat(sol.u...)
    # Perform convolution for each component
    ConvKin = Vector{Array{Float64,1}}(undef, num_states)
    for i in 1:num_states # loop through components
        component = sol_array[i,:]
        ConvKin[i] = GlobalAnalysis.convolveIRF(time_stamp, component, μ, σ, TimeIRFPos, IRFWindowCut, TimeODEBool)
    end

    KinMatrix = reduce(vcat, ConvKin')  # kinetic matrix
    
    TestSpec = data / KinMatrix  # spectral signatures matrix
    testData = TestSpec * KinMatrix # test matrix from kinetics and signatures

    return testData, TestSpec, KinMatrix
end

# objective function which returns the sum of squared differences bewteen 
# the produced data with guessed parameters (GetData) and the real data matrix
function Objective(param)
    da = GetData(param)
    testData = da[1]
    return sum(abs2, testData .- data)
end

# minimize objective function to optimize parameters
OP = bboptimize(Objective; SearchRange = bounds)
RecoveredParam = best_candidate(OP)

RecoveredData = GetData(RecoveredParam)
DataMatrix = RecoveredData[1] # data matrix
RecoveredSpec = RecoveredData[2] # spectral signatures
RecoveredKinetics = RecoveredData[3] # kinetic traces