# precompile workload for PackageCompiler, runs during sysimg building

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

cd(REPO_ROOT)
include(joinpath(REPO_ROOT, "src", "setup.jl"))

file = joinpath(REPO_ROOT, "data", "testData_first_order_seq.csv")
data = importData(file)[1]
s = data.x
t = data.y
d = data.z

t_bl = t .> 0
heatmap(t[t_bl], s, d[:, t_bl], xscale=:log10, c=:viridis,
    xlabel="Time (ps)", ylabel="Wavelength (nm)")

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

mh_bounds = generateBoundsMH(rn, limits)
_, fit_bounds, ode_helpers = setupVariables(rn, limits)
ssr_data = setupSSRMetaData(d)
param = vec((fit_bounds[:, 1] .+ fit_bounds[:, 2]) ./ 2)
batch = repeat(transpose(param), 4, 1)

paramToKin(t, param, ode_helpers)
paramToData(t, param, d, ode_helpers)
paramToSSR(t, param, d, ode_helpers, ssr_data)
paramToSSRParallel(t, batch, d, ode_helpers, ssr_data)

options = Options(iterations=1, parallel_evaluation=true, store_convergence=false)
Metaheuristics.optimize(
    p -> paramToSSRParallel(t, p, d, ode_helpers, ssr_data),
    mh_bounds,
    DE(N=4; options),
)
