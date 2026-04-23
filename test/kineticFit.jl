using Test
using Random
using LinearAlgebra
using Metaheuristics

include(joinpath(@__DIR__, "..", "src", "import.jl"))
include(joinpath(@__DIR__, "..", "src", "kinetic.jl"))

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

const sequentialFitReferenceKin = [
    0.023043777755121196 8.148492458492275e-5 2.7594832041848413e-6
    0.16023410936201843 0.0007955373525667625 3.5969236312177794e-5
    0.5064316548893468 0.003811772907254229 0.0002414005405283092
    0.9692364276184998 0.02504975016490156 0.0043442183814653545
    0.9227387492217802 0.05268157458436241 0.024507375130612005
    0.7561448745312697 0.07755566110601918 0.16439312484031207
    0.372931870293659 0.04168028543065749 0.5563703676716216
    0.05076514126213349 0.005671712716481053 0.7708358656866977
    4.1318669844827724e-5 4.618575714033562e-6 0.4113395884672498
    -9.114668339194349e-18 -1.0188778065957368e-18 0.05591643996260802
]

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
    return (; t, d, bounds, odeHelpers, ssrData)
end

@testset "Direct SSR Objective Matches Residual Path" begin
    problem = buildSequentialFitProblem()
    center = vec(mean(problem.bounds, dims=2))
    offcenter = center .* [1.05, 0.95, 1.02, 0.80, 1.10]

    for param in (center, offcenter)
        direct = paramToSSR(problem.t, param, problem.d, problem.odeHelpers, problem.ssrData)
        residuals = paramToResiduals(problem.t, param, problem.d, problem.odeHelpers)
        ssr = nansum(abs2.(residuals))
        @test isapprox(direct, ssr; rtol=1e-12, atol=1e-12)
    end
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

    return (; result, fitParam, fitKin)
end

@testset "Sequential Kinetic Fit Regression" begin
    regression = runSequentialFitRegression()
    fitParam = regression.fitParam
    sampledKin = regression.fitKin[sequentialFitSampleIdx, :]

    # objective stays in the good-fit basin
    @test minimum(regression.result) < 700.0
    # fitted parameters remain near the reference solution
    @test all(abs.(fitParam .- sequentialFitReferenceParam) .<= sequentialFitParamTolerance)
    # no sampled kinetic point drifts too far
    @test maximum(abs.(sampledKin .- sequentialFitReferenceKin)) < 0.1
    # overall sampled kinetics remain close in norm
    @test norm(sampledKin .- sequentialFitReferenceKin) / norm(sequentialFitReferenceKin) < 0.08
end
