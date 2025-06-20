
"""
    getNumFitVar(limits)

Return the number of fitted parameters.

# Arguments
- `limits` : dictionary containing starting values and bounds for 
the parameter optimization   

"""
function getNumFitVar(limits)
    # count fit parameters
    count = 0
    for (k,v) in limits
        # if dict entry is a 2-element vector or tuple, it is a
        # fitted variable (fixed variables have a single value)
        if length(v) == 2
            count += 1
        end
    end

    return count
end


"""
    calculateAdjR2(yObs, yFit, k)

Return the adjusted R² for a model with `k` fitted parameters.

# Arguments
- `yObs` : array of observed values  
- `yFit` : array of model-predicted values
- `k`    : number of free parameters in the fitted model

"""
function calculateAdjR2(yObs, yFit, k)

    # flatten matrices to vectors
    yObs = vec(yObs)
    yFit = vec(yFit)

    # experimental and fit vectors much have same length
    @assert length(yObs) == length(yFit)

    n = length(yObs)
    # sum of squared residuals
    ssr = nansum((yObs .- yFit) .^ 2)
    ȳ = nanmean(yObs)
    # total sum of squares
    sst = nansum((yObs .- ȳ) .^ 2)

    R² = 1 - ssr / sst
    adjR² = 1 - (1 - R²) * (n - 1) / (n - k - 1)

    return adjR²
end


"""
    chiSq(yObs, yFit; σ = 1, k = 0)

Compute the Pearson chi-square statistic and its reduced form.

Returns a tuple `(χ², redχ²)` where  

* `χ²` : chi-squared
* `redχ²` : reduced chi-squared

# Arguments
- `yObs` : array of observed values  
- `yFit` : array of model-predicted values
- `k` : number of free parameters in the model (used only for `redχ²`).
- `σ` : per-point standard deviation(s). A scalar is broadcast to all 
points; a vector must match `yObs` length. Default `σ = 1` gives the 
unweighted χ².  

"""
function calculateChiSq(yObs, yFit, k; σ=1)

    # flatten matrices to vectors
    yObs = vec(yObs)
    yFit = vec(yFit)

    # experimental and fit vectors much have same length
    @assert length(yObs) == length(yFit)

    σVec = isa(σ, AbstractArray) ? vec(σ) : fill(σ, length(yObs))
    # all σ values must be positive
    @assert all(σVec .> 0)

    # calculate chi-squared
    residuals = (yObs .- yFit) ./ σVec
    χ² = nansum(residuals .^ 2)

    # calculate reduced chi-squared
    dof = length(yObs) - k
    redχ² = χ² / dof

    return χ², redχ²
end


"""
    calculateAICc(yObs, yFit, k)

Return AICc for a least-squares fit.

# Arguments
- `yObs` : array of observed values  
- `yFit` : array of model-predicted values
- `k` : number of free parameters in the model

AICc is reliable when the residuals are independent, Gaussian-like, 
and the same data set is used across models; it remains unbiased for 
small samples, but loses meaning if parameters are unidentifiable or 
the noise model is misspecified.

"""
function calculateAICc(yObs, yFit, k)
    
    # flatten matrices to vectors
    yObs = vec(yObs)
    yFit = vec(yFit)

    # experimental and fit vectors much have same length
    @assert length(yObs) == length(yFit)

    n = length(yObs)
    # sum of squared residuals
    ssr = nansum((yObs .- yFit) .^ 2)
    aic = 2k + n*log(ssr/n)
    aicc = aic + (2k*(k+1)) / (n - k - 1)

    return aicc
end


"""
    calculateBIC(yObs, yFit, k)

Return the Bayesian Information Criterion (BIC) for a least-squares fit.

# Arguments
- `yObs` : array of observed values  
- `yFit` : array of model-predicted values
- `k` : number of free parameters in the model

BIC assumes the candidate set contains the true model and penalises complexity 
more strongly; its χ² + k ln n formula is asymptotically valid only for large n, 
regular (well-identified) parameters, and independent Gaussian errors.

"""
function calculateBIC(yObs, yFit, k)
    # flatten matrices to vectors
    yObs = vec(yObs)
    yFit = vec(yFit)

    # experimental and fit vectors much have same length
    @assert length(yObs) == length(yFit)

    n = length(yObs)
    # sum of squared residuals
    ssr = nansum((yObs .- yFit) .^ 2)
    bic = k*log(n) + n*log(ssr/n)

    return bic
end
