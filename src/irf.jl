using Distributions
using Interpolations 
using DSP
using MonteCarloMeasurements


#----------------------------------------------------------------------
# Generate optimised time vectors to calculate kinetic traces over
#----------------------------------------------------------------------

"""
Returns unique steps in vector `t`.
"""
function getTimeSteps(t)
    # round to avoid machine precision inaccuracies
    unique(round.(diff(t); digits=14))
end


"""
Replace particle-valued IRF parameters by their mean values.

The current IRF implementation uses a deterministic time grid and Gaussian
kernel, so particle information in `μ` and `σ` is intentionally collapsed
before grid construction and convolution. This approach will be reevaluated
later.
"""
function getDeterministicIrfParameters(μ, σ)
    if μ isa AbstractParticles
        μ = pmean(μ)
    end
    if σ isa AbstractParticles
        σ = pmean(σ)
    end
    return μ, σ
end


"""
Returns whether `t` has a constant step, the first rounded step, and the
minimum rounded step, while avoiding allocation of the full `diff(t)`.
"""
function getTimeStepInfo(t)
    @assert length(t) ≥ 2

    firstStep = round(t[2] - t[1]; digits=14)
    minStep = firstStep
    isConstant = true

    @inbounds for n in 3:length(t)
        step = round(t[n] - t[n-1]; digits=14)
        if step != firstStep
            isConstant = false
        end
        if step < minStep
            minStep = step
        end
    end

    return isConstant, firstStep, minStep
end


"""
Merge sorted vectors `a` and `b` into a sorted unique vector and return the
positions of the original entries in that merged vector.
"""
function mergeSortedUniqueWithIndices(a, b)
    nA = length(a)
    nB = length(b)
    merged = Vector{Float64}(undef, nA + nB)
    aIdx = Vector{Int}(undef, nA)
    bIdx = Vector{Int}(undef, nB)

    i = 1
    j = 1
    k = 1
    while i ≤ nA && j ≤ nB
        ai = a[i]
        bj = b[j]
        if isapprox(ai, bj; atol=1e-12, rtol=1e-12)
            merged[k] = ai
            aIdx[i] = k
            bIdx[j] = k
            i += 1
            j += 1
            k += 1
        elseif ai < bj
            merged[k] = ai
            aIdx[i] = k
            i += 1
            k += 1
        else
            merged[k] = bj
            bIdx[j] = k
            j += 1
            k += 1
        end
    end

    while i ≤ nA
        merged[k] = a[i]
        aIdx[i] = k
        i += 1
        k += 1
    end

    while j ≤ nB
        merged[k] = b[j]
        bIdx[j] = k
        j += 1
        k += 1
    end

    resize!(merged, k - 1)
    return merged, aIdx, bIdx
end


"""
Build a mirrored grid `[-reverse(tPos); tPos]` without temporary arrays.
"""
function buildMirroredTimeGrid(tPos)
    n = length(tPos)
    out = Vector{Float64}(undef, 2n)
    @inbounds for i in 1:n
        out[i] = -tPos[n - i + 1]
        out[n + i] = tPos[i]
    end
    return out
end


"""
Returns time vector for ODE solver `tOde` given an input `t`
in which all time points have equal spacing. 

Extends `t` by additional points, which will be discarded after 
calculating the kinetic trace to avoid edge effect over the 
original time range.
"""
function getOdeTimeConstantStep(t, μ, σ)
    _, tStep, _ = getTimeStepInfo(t)

    # final extended time point
    tExtFinal = t[end]+4σ-μ
    # index of first positive point in t
    firstPos = searchsortedlast(t, 0.0) + 1
    nBase = max(0, length(t) - firstPos + 1)
    tExt = (t[end] + tStep):tStep:tExtFinal
    nExt = length(tExt)

    tOde = Vector{Float64}(undef, nBase + nExt)
    if nBase > 0
        copyto!(tOde, 1, t, firstPos, nBase)
    end
    @inbounds for n in 1:nExt
        tOde[nBase + n] = tExt[n]
    end

    return (tOde,)

end


"""
Returns time vector for ODE solver `tOde` given an input `t`
in which time points have arbitrary spacing.

`tOde` contains early time points with equal spacing, used for convolution
with the instrument response, and arbitrarily spaced points at later
times where effect of instrument response is negligible.
"""
function getOdeTimeVariableStep(t, μ, σ)
    isConstant, _, minStep = getTimeStepInfo(t)
    @assert !isConstant

    # time step in irf window is smallest time step in dataset
    irfStep = minStep
    # IRF window width in pos/neg direction
    tMaxIrf = 30*σ
    # evenly spaced time steps for convolution
    # need one point after tMaxIrf in case interpolation is queried between final step and tMaxIrf
    tConv = getOdeTimeConstantStep(irfStep/2:irfStep:tMaxIrf+irfStep, μ, σ)[1]

    # time steps post convolution; shift post-irf time by -μ to 
    # evaluate correct points for original time vector in ODE
    firstPost = searchsortedfirst(t, tMaxIrf)
    nPost = max(0, length(t) - firstPost + 1)
    tPostConv = Vector{Float64}(undef, nPost)
    @inbounds for n in 1:nPost
        tPostConv[n] = t[firstPost + n - 1] - μ
    end
    
    # combine irf and post-irf times in array for single call to ODE solver
    tOde, tConvIdx, tPostIdx = mergeSortedUniqueWithIndices(tConv, tPostConv)
    
    return (tOde, tConv, tMaxIrf, tConvIdx, tPostIdx)

end


"""
Returns vector over which to calculate kinetic traces
based on `t`, optimised for a Gaussian IRF with centre `μ` and
standard deviation `σ`, and automatically taking care of 
`t` spacing. 
"""
function getOdeTime(t, μ, σ)
    μ, σ = getDeterministicIrfParameters(μ, σ)

    isConstant, _, _ = getTimeStepInfo(t)

    # constant time spacing
    if isConstant
        return getOdeTimeConstantStep(t, μ, σ)
    # variable time spacing
    else
        return getOdeTimeVariableStep(t, μ, σ)
    end
end


#----------------------------------------------------------------------
# Convolution of kinetic trace(s) and IRF
#----------------------------------------------------------------------


"""
Performs numerical integration of `y` over `x` via trapezoidal method.
"""
function trapezIntegration(x, y) 
    # x and y must have same length
    @assert length(x) == length(y)
    out = 0.0
    for n in 2:length(x)
      out += 0.5*(x[n] - x[n-1])*(y[n] + y[n-1])
    end
    return out
end


"""
Generates Gaussian over vector `t`, centred around `μ` with 
standard deviation `σ`.
"""
function getGaussianIRF(t, μ, σ)
    μ, σ = getDeterministicIrfParameters(μ, σ)
    irf = Vector{Float64}(undef, length(t))
    norm = inv(σ * √(2π))
    @inbounds for n in eachindex(t)
        d = (t[n] - μ) / σ
        irf[n] = norm * exp(-0.5 * d * d)
    end
    return irf
end


"""
Convolves vectors `kin` and `irf`, where `kin` contains 
data for `t ≥ 0` only and must be evenly spaced in `t`. 
`t` contains positive and negative time points, `irf` 
is calculated over full `t`.

`kin` can be a single kinetic trace or a matrix consisting of
multiple kinetic traces as columns.
"""
function convolveIrfConstantStep(t, kin, irf, tStep::Real)
    @assert length(t) == length(irf)

    # normalise IRF by its area
    irf ./= trapezIntegration(t,irf)

    # discrete convolution, hence multiply by time step
    kinConv = DSP.conv(kin, irf)
    kinConv .*= tStep

    # make convolved data same length as IRF data to restore 
    # time correspondence
    kinConvSame = selectdim(kinConv, 1, 1:length(t))

    return kinConvSame

end

"""
Compatibilty function in case `tSteps` is an array.
"""
function convolveIrfConstantStep(t, kin, irf, tSteps) 
    @assert length(tSteps) == 1
    return convolveIrfConstantStep(t, kin, irf, tSteps[1])
end



"""
Convolves vectors `kin` and a Gaussian IRF with centre `μ` and
standard deviation `σ`. `kin` is contains data for `t ≥ 0` only 
and can be arbitrarily spaced in `t`. `t` contains positive and 
negative time points. `tStepParam` is generated by `getOdeTime`.

`kin` can be a single kinetic trace or a matrix consisting of
multiple kinetic traces as columns. These kinetic traces can be
of arbitrary reaction order.

To convolve a kinetic traces with variable time spacing, the kinetic
trace is calculated over an evenly spaced early time period and then
interpolated back onto the early portion of the original `t` vector. 
Later times, at which the effect of the instrument response is 
negiligible, are appended to the convolved data without convolution.
"""
function convolveIrfVariableStep(t, kin, μ, σ, tSteps, tStepParam)
    tOde, tConv, tMaxIrf, tConvIdx, tPostIdx = tStepParam

    # mirror evenly spaced times around zero
    tIrf = buildMirroredTimeGrid(tConv)
    # generate IRF
    irf = getGaussianIRF(tIrf, μ, σ)

    # convolve evenly spaced part of kinetic trace
    irfStep = length(tSteps) == 1 ? tSteps[1] : tConv[2] - tConv[1]
    if ndims(kin) == 1
        kinConv = convolveIrfConstantStep(tIrf, view(kin, tConvIdx), irf, irfStep)
    else
        kinConv = convolveIrfConstantStep(tIrf, view(kin, tConvIdx, :), irf, irfStep)
    end

    activeStart = searchsortedfirst(t, -tMaxIrf)
    postStart = searchsortedfirst(t, tMaxIrf)

    # interpolator for convolved trace
    # single kinetic trace
    if ndims(kin) == 1
        itpKinConv = interpolate((tIrf,), kinConv, Gridded(Linear()))
        kinConvFinal = similar(kin, length(t))
        fill!(kinConvFinal, zero(eltype(kinConvFinal)))
        @inbounds for i in activeStart:postStart-1
            kinConvFinal[i] = itpKinConv(t[i])
        end
        @inbounds for i in postStart:length(t)
            kinConvFinal[i] = kin[tPostIdx[i - postStart + 1]]
        end
    # multiple kinetic traces
    else
        # no interpolation along second dimension 
        numKin = size(kin,2)
        itpKinConv = interpolate((tIrf,1:numKin,), kinConv, 
            (Gridded(Linear()),NoInterp()))
        kinConvFinal = similar(kin, length(t), numKin)
        fill!(kinConvFinal, zero(eltype(kinConvFinal)))
        @inbounds for m in 1:numKin
            for i in activeStart:postStart-1
                kinConvFinal[i,m] = itpKinConv(t[i], m)
            end
            for i in postStart:length(t)
                kinConvFinal[i,m] = kin[tPostIdx[i - postStart + 1],m]
            end
        end
    end

    return kinConvFinal

end


"""
Convolves vectors `kin` and a Gaussian IRF with centre `μ` and
standard deviation `σ`, automatically taking care of `t` spacing.
For variable `t` spacing, `tStepParam` must be supplied from
`getOdeTime`.
"""
function convolveIRF(t, kin, μ, σ, tStepParam)

    isConstant, tStep, minStep = getTimeStepInfo(t)

    # constant time spacing
    if isConstant    
        # generate IRF
        irf = getGaussianIRF(t, μ, σ)
        return convolveIrfConstantStep(t, kin, irf, tStep) 
    # variable time spacing
    else
        return convolveIrfVariableStep(t, kin, μ, σ, (minStep,), tStepParam)
    end

end



