using Test
using Random
using LinearAlgebra
using Statistics
using Metaheuristics
using Catalyst
using GlobalAnalysis
using NaNStatistics

const sequentialFitReferenceParam = [
    0.9977211681177615,
    0.1003118642822456,
    0.010003001737107244,
    0.19957692736217422,
    0.0997021182222694,
]

const sequentialFitParamTolerance = [
    0.15,
    0.02,
    0.004,
    0.03,
    0.04,
]

const sequentialFitSampleIdx = [101, 102, 143, 196, 237, 301, 372, 436, 507, 571]

function buildSequentialFitProblem()
    file = joinpath(@__DIR__, "..", "data", "testData_first_order_seq.csv")
    data = importData(file)[1]
    t = data.y
    d = data.z

    rn = @reaction_network begin
        k1, A --> B
        k2, B --> C
        k3, C --> 0
    end

    limits = Dict(
        :A => 1,
        :B => 0,
        :C => 0,
        :k1 => [5e-1, 5],
        :k2 => [5e-2, 5e-1],
        :k3 => [5e-3, 5e-2],
        :μ => [-0.5, 0.5],
        :σ => [0.04, 0.2],
    )

    _, bounds, odeHelpers = setupVariables(rn, limits)
    ssrData = setupSSRMetaData(d)
    species = odeHelpers[5]
    return (; t, d, bounds, odeHelpers, ssrData, species)
end

@testset "Direct SSR Objective Matches Residual Path" begin
    problem = buildSequentialFitProblem()
    center = vec(mean(problem.bounds, dims=2))
    offcenter = center .* [1.05, 0.95, 1.02, 0.80, 1.10]

    for param in (center, offcenter)
        direct = paramToSSR(problem.t, param, problem.d, problem.odeHelpers, problem.ssrData)
        residuals = paramToResiduals(problem.t, param, problem.d, problem.odeHelpers)
        ssr = nansum(abs2.(residuals))
        @test isapprox(direct, ssr; rtol=1e-6, atol=1e-6)
    end
end

@testset "Parallel SSR Objective Handles Changing IRF Grids" begin
    problem = buildSequentialFitProblem()
    lo = problem.bounds[:, 1]
    hi = problem.bounds[:, 2]
    batch = Matrix{Float64}(undef, 24, length(lo))

    for n in axes(batch, 1)
        frac = (n - 1) / (size(batch, 1) - 1)
        batch[n, :] .= lo .+ frac .* (hi .- lo)
    end

    ssr = paramToSSRParallel(problem.t, batch, problem.d, problem.odeHelpers, problem.ssrData)

    @test length(ssr) == size(batch, 1)
    @test all(isfinite, ssr)
end

function runSequentialFitRegression(; seed=1234, iterations=30, population=8)
    problem = buildSequentialFitProblem()
    objective = param -> paramToSSR(problem.t, param, problem.d, problem.odeHelpers, problem.ssrData)

    blasThreads = BLAS.get_num_threads()
    result = nothing
    try
        BLAS.set_num_threads(1)

        Random.seed!(seed)
        warmupOptions = Options(iterations=2, store_convergence=false)
        Metaheuristics.optimize(objective, problem.bounds, DE(N=population; options=warmupOptions))

        Random.seed!(seed)
        options = Options(iterations=iterations, store_convergence=false)
        result = Metaheuristics.optimize(objective, problem.bounds, DE(N=population; options=options))
    finally
        BLAS.set_num_threads(blasThreads)
    end

    fitParam = minimizer(result)
    _, _, fitKin = paramToData(problem.t, fitParam, problem.d, problem.odeHelpers)
    _, _, referenceKin = paramToData(problem.t, sequentialFitReferenceParam, problem.d, problem.odeHelpers)

    return (; result, fitParam, fitKin, referenceKin, species=problem.species)
end

@testset "Sequential Kinetic Fit Regression" begin
    regression = runSequentialFitRegression()
    fitParam = regression.fitParam
    referenceSpecies = [:A, :B, :C]
    speciesOrder = [findfirst(==(species), regression.species) for species in referenceSpecies]
    @test all(!isnothing, speciesOrder)
    sampledKin = regression.fitKin[sequentialFitSampleIdx, Int.(speciesOrder)]
    sampledReferenceKin = regression.referenceKin[sequentialFitSampleIdx, Int.(speciesOrder)]

    # objective stays in the good-fit basin
    @test minimum(regression.result) < 700.0
    # fitted parameters remain near the reference solution
    @test all(abs.(fitParam .- sequentialFitReferenceParam) .<= sequentialFitParamTolerance)
    # no sampled kinetic point drifts too far
    @test maximum(abs.(sampledKin .- sampledReferenceKin)) < 0.1
    # overall sampled kinetics remain close in norm
    @test norm(sampledKin .- sampledReferenceKin) / norm(sampledReferenceKin) < 0.08
end
