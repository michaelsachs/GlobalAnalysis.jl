using Catalyst
using DifferentialEquations
using NaNStatistics
using Distributions
using FiniteDiff
using LinearAlgebra
using MonteCarloMeasurements

include("irf.jl")


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
    p = Dict{Symbol,eltype(fitParam)}()
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
    species = Catalyst.unknowns(rn)
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
function paramToData(t, rn, param, limits, Data)
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

    # switch solver depending on data type
    if eltype(param) isa Float64
        sol  = solve(prob, AutoTsit5(Rosenbrock23()))
    else
        # for MonteCarloMeasurements
        sol  = solve(prob,  Rosenbrock23(autodiff=AutoFiniteDiff()))
    end

    kin = transpose(Array(sol))

    # convolve kinetic traces with IRF
    kinConv = convolveIRF(t, kin, μ, σ, tStepParam)

    # generate spectra based on calculated kinetics
    testSpc = Data / kinConv'
    # assemble data matrix 
    testData = testSpc * kinConv'

    return testData, testSpc, kinConv

end


"""
Returns 2D matrix of residuals, calculated by subtracting the simulated 
matrix from the experimental one.
"""
function paramToResiduals(t, rn, param, limits, Data)
    testData, _, _ = paramToData(t, rn, param, limits, Data)
    return testData .- Data
end


"""
Returns vector of residuals by flattening the output of `paramToResiduals`.
Used for jacobian calculation.
"""
function paramToResidualsVec(t, rn, param, limits, Data)
    res = paramToResiduals(t, rn, param, limits, Data)
    return vec(res)
end


"""
Returns sum of squared residuals between simulated and experimental
data. Used for parameter optimization.
"""
function paramToSSR(t, rn, param, limits, Data)
    res = paramToResiduals(t, rn, param, limits, Data)
    return nansum((res).^2)
end


"""
Parallel evaluation of `paramToSSR` for use with Metaheuristics.
`param` is a 2D array, with dimension 1 corresponding to the number
of parallel evaluations, and dimension 2 being the number of fit
parameters.
"""
function paramToSSRParallel(t, rn, param, limits, Data)
    ssr = zeros(size(param,1))
    Threads.@threads for n in 1:size(param,1)
        ssr[n] = paramToSSR(t, rn, param[n,:], limits, Data)
    end
    return ssr
end


"""
Generates string array to label species in reaction network `rn`
in plots.
"""
function getLabels(rn)
    permutedims(string.(getSpecies(rn)))
end


"""
Returns the half-widths `halfCI` of the two-sided confidence intervals for the
best-fit parameter vector `paramOpt`.  The interval for each parameter is

    paramOpt[i] ± halfCI[i]

and is expected to contain the true value in roughly `confidenceLevel*100%`  
of repeated experiments (e.g. `confidenceLevel=0.95` → 95% confidence interval).

Uncertainties are computed from the QR-based Gaussian approximation.
"""
function getParamConfidence(t, rn, paramOpt, limits, Data; confidenceLevel=0.95)

    resVec = paramToResidualsVec(t, rn, paramOpt, limits, Data)
    # number of observation (time x energy)
    numObs = length(resVec)
    # degrees of freedom
    dof = numObs-length(paramOpt)
    # sum of squared residuals
    ssr = nansum(resVec.^2)

    # make sure that we have more observables than parameters (ensure 
    # tall-matrix Jacobian, without which R is not invertible)
    if dof ≤ 0
        error("The number of parameters exceeds the number of observables.")
    end

    # calculate the jacobian matrix
    J = FiniteDiff.finite_difference_jacobian(
            p->paramToResidualsVec(t, rn, p, limits, Data), paramOpt)
    # perform QR decomposition
    F = LinearAlgebra.qr(J)
    # get upper triangular matrix R
    R = F.R
    # calculate R-inverse
    Rinv = LinearAlgebra.inv(R)

    # calculate noise-variance estimate
    σ² = ssr / dof
    # calculate full variance–covariance matrix
    cov = σ² * (Rinv * Rinv')
    # get full parameter variances
    v = diag(cov)

    # standard error for each parameter; one-sigma (68%) 
    stError = sqrt.(v)

    # two-sided tail area (e.g. 0.025 for 95%)
    α = (1-confidenceLevel)/2
    # Student-t critical value for those tails
    tScore = quantile(TDist(dof), 1 - α)
    # half-width of the (1–α) confidence interval
    halfCI = tScore .* stError

    return halfCI, cov

end


"""
Calculates confidence intervals on kinetics and spectra obtained from optimized parameters
using Monte Carlo sampling. `samples` is the number of Monte Carlo samples and 
`confidenceLevel` is the desired confidence level.
"""
function paramToDataCI(t, rn, paramOpt, limits, Data; samples=200, confidenceLevel=0.95)

    # get parameter covariance matrix at the desired confidence level
    _,cov = getParamConfidence(t, rn, paramOpt, limits, Data; confidenceLevel=confidenceLevel)

    # create a multivariate normal distribution
    sampler = MvNormal(paramOpt, cov)
    # vector of particle objects
    particles = StaticParticles(samples, sampler)

    # propagate particles through the forward model
    fitData,fitSpc,fitKin = paramToData(t, rn, particles, limits, Data)

    return fitData,fitSpc,fitKin

end





