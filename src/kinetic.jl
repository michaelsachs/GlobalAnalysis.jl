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
    setupVariables(rn, limits)

Pre-process the reaction network once and return everything the
optimiser and forward model need.

Inputs
------
* `rn`: Catalyst reaction network  
* `limits`  
  * `[:x] => v`         # fixed at value v
  * `[:x] => [lo, hi]`  # fitted, bounded lo-hi

Outputs
-------
* `syms` : all variables in a fixed order: `[species; rateConstants; μ; σ]`
* `fitBounds` : `n×2` array of lower/upper bounds for the free parameters,
  in the exact order the optimiser will see them.
* `odeHelpers` : helper vector reused by the forward model:

  1. `paramTempl` — complete parameter vector with fixed variables
     inserted and zeros for fitted variables.
  2.  3-element vector of vectors containing indices of fit variables:  
     `fitIdx[1]` fitted-species indices,  
     `fitIdx[2]` fitted-rate-constant indices,  
     `fitIdx[3]` fitted-IRF indices.
  3. `fitLinearIdx` — flattened version of `fitIdx`.
  4. `idxRanges`  — the three contiguous index ranges that correspond to
     *all* species, rate constants, and IRF parameters respectively.
  5. `species`    — species symbols, used when particle-valued ODE inputs
     need dictionary-based initial conditions.
  6. `rateConst`  — rate-constant symbols, used when particle-valued ODE
     inputs need dictionary-based parameters.
  7. `rn`         — reaction network used to build the ODE problem.

Notes
-----
* All dictionary look-ups and symbol handling are done here, once,
  so the threaded objective works mostly with numeric vectors and
  precomputed helper metadata.
"""
function setupVariables(rn, limits)

    species = getSpecies(rn)
    rateConst = getParameters(rn)
    irfParam = [:μ; :σ]

    syms = [species; rateConst; irfParam]
    # parameter template with values for fixed variables and 
    # zeroes for fit variables
    paramTempl = zeros(Float64, length(syms))
    # true indicates that variable is fitted
    fitMask = falses(length(syms))
    # gather fit bounds
    fitBounds = []

    for (n,sym) in pairs(syms)
        if haskey(limits, sym)
            val = limits[sym]
            # fixed variable
            if length(val) == 1
                paramTempl[n] = val
            # fit variable
            elseif length(val) == 2
                fitMask[n] = true
                push!(fitBounds, [val...])
            else
                error("Bounds value for $(sym) must be
                (a) a number if $(sym) is a fixed parameter, or
                (b) a vector containing two numbers if $(sym) is a fit parameter")
            end
        else
            error("Variable $(sym) is part of reaction network but missing in limits dictionary")
        end
    end

    # transforms bounds into nx2 array so each row matches one fit parameter
    fitBounds = permutedims(hcat(fitBounds...))
    # indices for fit variables
    fitIdx = findall(fitMask)

    lengths = length.([species, rateConst, irfParam])
    # starting indices
    firsts = cumsum(vcat(1, lengths[1:end-1]))
    # index ranges for each component
    idxRanges = [s:(s+l-1) for (s,l) in zip(firsts, lengths)]
    # indices of fitted variables for each component 
    fitIdx = [ [i for i in grp if fitMask[i]]  for grp in idxRanges ]
    # flatten indices
    fitLinearIdx = vcat(fitIdx...)

    odeHelpers = [paramTempl, fitIdx, fitLinearIdx, idxRanges, species, rateConst, rn]

    return syms, fitBounds, odeHelpers

end


"""
Generates convolved kinetic traces based on time vector `t`, reaction
network `rn`, and parameter vector `param`. `param` contains parameters
in order (1) initial state populations, (2) rate constants, (3) IRF
parameters.
"""
function paramToKin(t, param, odeHelpers)

    # split helper array into components
    paramTempl, _, fitLinearIdx, idxRanges, species, rateConst, rn = odeHelpers

    pType = eltype(param)

    if pType == Float64
        _paramTempl = copy(paramTempl)
    else
        # 1 ± 0 as Particles
        base = one(pType)          
        _paramTempl = [paramTempl[i] * base for i in eachindex(paramTempl)]
    end

    # fill fit parameters into template
    @inbounds _paramTempl[fitLinearIdx] .= param
    # distribute parameters from template
    u0, ks, irf = view.(Ref(_paramTempl), idxRanges)

    # convert vector to dict; essential for ODE solver to work with Particles
    if pType ≠ Float64
        u0 = Dict(species[i] => u0[i] for i in eachindex(u0))
        ks = Dict(rateConst[i] => ks[i] for i in eachindex(ks))
    end

    # assemble time vector for ODE solver
    tStepParam = getOdeTime(t, irf...)
    tOde = tStepParam[1]
    # time span for ODE solver
    tspan = [minimum(tOde), maximum(tOde)]

    # set up and solve ODEs
    prob = ODEProblem(rn, u0, tspan, ks; saveat=tOde, save_everystep=false, dense=false)

    # switch solver depending on data type
    if pType == Float64
        sol  = solve(prob, AutoTsit5(Rosenbrock23()))
    # for MonteCarloMeasurements 
    else
        # Rosenbrock23 requires unsafe comparisons for Particles, avoid for now
        sol = solve(prob, Tsit5())
    end

    kin = transpose(Array(sol))

    # convolve kinetic traces with IRF
    return convolveIRF(t, kin, irf..., tStepParam)

end


"""
Generates test data by calculating kinetic traces based on time vector
`t`, reaction network `rn`, and parameter vector `param`. `param` contains 
parameters in order (1) initial state populations, (2) rate constants, 
(3) IRF parameters.
"""
function paramToData(t, param, Data, odeHelpers)

    kinConv = paramToKin(t, param, odeHelpers)

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
function paramToResiduals(t, param, Data, odeHelpers)
    testData, _, _ = paramToData(t, param, Data, odeHelpers)
    return testData .- Data
end


"""
Returns vector of residuals by flattening the output of `paramToResiduals`.
Used for jacobian calculation.
"""
function paramToResidualsVec(t, param, Data, odeHelpers)
    res = paramToResiduals(t, param, Data, odeHelpers)
    return vec(res)
end


"""
Precompute SSR metadata for a fixed dataset so repeated objective evaluations
avoid rescanning `Data`.
"""
function setupSSRMetaData(Data)
    dataNorm = zero(eltype(Data))
    hasNaN = false

    @inbounds for datum in Data
        if isnan(datum)
            hasNaN = true
        else
            dataNorm += datum * datum
        end
    end

    return SSRMetaData(dataNorm, hasNaN)
end


"""
Fastest SSR path which avoids constructing residuals column-by-column. Requires
input data without NaNs.

`testData = testSpc * kinConv'` is the closest reconstruction of `Data`
reachable with the current basis (`kinConv`). The remaining residual
`Data - testData` is therefore pure unexplained error, so:
``||Data - testData||² = ||Data||² - ||testData||²``.
"""
function projectSSR(dataNorm, testSpc, kinConv)
    gram = kinConv' * kinConv
    projectedSpc = testSpc * gram
    ssr = dataNorm - dot(testSpc, projectedSpc)

    # Guard against tiny negative values from floating-point roundoff.
    if ssr < 0 && isapprox(ssr, zero(ssr); atol=eps(real(ssr)) * max(one(real(ssr)), dataNorm))
        return zero(ssr)
    end

    return ssr
end


"""
NaN-tolerant in-place residual accumulation that reuses the column buffer
`_spc`.
"""
function accumulateSSR!(_spc, Data, testSpc, kinConv)

    @assert size(Data, 1) == size(testSpc, 1)
    @assert size(Data, 2) == size(kinConv, 1)
    @assert size(testSpc, 2) == size(kinConv, 2)

    # create a correctly-typed zero for accumulating ssr
    ssr = zero(promote_type(eltype(Data), eltype(testSpc), eltype(kinConv)))

    @inbounds for n in axes(Data, 2)
        mul!(_spc, testSpc, view(kinConv, n, :))
        for m in axes(Data, 1)
            diff = _spc[m] - Data[m, n]
            if !isnan(diff)
                ssr += diff * diff
            end
        end
    end

    return ssr

end


"""
Returns sum of squared residuals between simulated and experimental
data. Used for parameter optimization.

`ssrData` must be precomputed once for the fixed dataset via and 
reused across objective evaluations.
"""
function paramToSSR(t, param, Data, odeHelpers, ssrData::SSRMetaData)
    kinConv = paramToKin(t, param, odeHelpers)
    testSpc = Data / kinConv'

    if ssrData.hasNaN
        # slower path if Data contains NaN
        _spc = similar(testSpc, size(Data, 1))
        ssr = accumulateSSR!(_spc, Data, testSpc, kinConv)
    else
        # fast path for Data without NaN
        ssr = projectSSR(ssrData.dataNorm, testSpc, kinConv)
    end

    return ssr
end


"""
Parallel evaluation of `paramToSSR` for use with Metaheuristics.
`param` is a 2D array, with dimension 1 corresponding to the number
of parallel evaluations (batch size), and dimension 2 being the number 
of fit parameters. This function is called once per fit iteration.
"""
function paramToSSRParallel(t, param, Data, odeHelpers, ssrData::SSRMetaData)
    ssr = zeros(size(param,1))
    Threads.@threads for n in 1:size(param,1)
        ssr[n] = paramToSSR(t, param[n,:], Data, odeHelpers, ssrData)
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
best-fit parameter vector `param`.  The interval for each parameter is

    paramOpt[i] ± halfCI[i]

and is expected to contain the true value in roughly `confidenceLevel*100%`  
of repeated experiments (e.g. `confidenceLevel=0.95` → 95% confidence interval).

Uncertainties are computed from the QR-based Gaussian approximation.
"""
function getParamConfidence(t, param, Data, odeHelpers; confidenceLevel=0.95)

    resVec = paramToResidualsVec(t, param, Data, odeHelpers)
    # number of observation (time x energy)
    numObs = length(resVec)
    # degrees of freedom
    dof = numObs-length(param)
    # sum of squared residuals
    ssr = nansum(resVec.^2)

    # make sure that we have more observables than parameters (ensure 
    # tall-matrix Jacobian, without which R is not invertible)
    if dof ≤ 0
        error("The number of parameters exceeds the number of observables.")
    end

    # calculate the jacobian matrix
    J = FiniteDiff.finite_difference_jacobian(
            p->paramToResidualsVec(t, p, Data, odeHelpers), param)
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
function paramToDataCI(t, param, Data, odeHelpers; samples=200, confidenceLevel=0.95)

    # get parameter covariance matrix at the desired confidence level
    _,cov = getParamConfidence(t, param, Data, odeHelpers; confidenceLevel=confidenceLevel)

    # create a multivariate normal distribution
    sampler = MvNormal(param, cov)
    # vector of particle objects
    particles = StaticParticles(samples, sampler)

    # propagate particles through the forward model
    fitData,fitSpc,fitKin = paramToData(t, particles, Data, odeHelpers)

    return fitData,fitSpc,fitKin

end


"""
Prints the optimized value of each fitted variable as its mean value `fitParam` ± confidence 
interval `ci`.
"""
function printFitResult(fitParam, ci, syms, fitLinearIdx)

    # vector of fit variables
    fitVar = syms[fitLinearIdx]
    # print each fit variable with ci
    for n in eachindex(fitVar)
        println("$(fitVar[n]): $(fitParam[n]) ± $(ci[n])")
    end

end
