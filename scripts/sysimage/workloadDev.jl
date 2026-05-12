# precompile workload for PackageCompiler, runs during sysimg building

using CSV
using Catalyst
using CommonSolve: solve
using DSP
using DelimitedFiles
using Distributions
using Downloads
using FiniteDiff
using Glob
using IJulia
using Interpolations
using LinearAlgebra
using Metaheuristics
using MonteCarloMeasurements
using NLopt
using NaNStatistics
using OrdinaryDiffEqRosenbrock: Rosenbrock23
using OrdinaryDiffEqTsit5: AutoTsit5
using Plots
using SciMLBase
using SciMLStructures
using Statistics
using StructArrays
using TOML
using Tables

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

cd(REPO_ROOT)

file = joinpath(REPO_ROOT, "data", "testData_first_order_seq.csv")
tbl = CSV.File(file)
mat = Tables.matrix(tbl)

heatmap(mat[2:5, 1:5], c=:viridis)

rn = @reaction_network begin
    k1, A --> B
    k2, B --> C
    k3, C --> 0
end

tspan = (0.0, 10.0)
u0 = [:A => 1.0, :B => 0.0, :C => 0.0]
p = [:k1 => 1.0, :k2 => 0.2, :k3 => 0.02]
prob = ODEProblem(rn, u0, tspan, p)
solve(prob, AutoTsit5(Rosenbrock23()); saveat=range(0.0, 10.0; length=16))

objective = x -> sum(abs2, x)
bounds = [-1.0 -1.0; 1.0 1.0]
options = Options(iterations=1, store_convergence=false)
Metaheuristics.optimize(objective, bounds, DE(N=4; options))
