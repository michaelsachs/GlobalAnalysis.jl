using PackageCompiler
using TOML

include("common.jl")

const HOST_CPU_FEATURES_PREFS = "HostCPUFeatures"
const HOST_CPU_FEATURES_UUID = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"

const DEV_PACKAGES = [
    :CSV,
    :Catalyst,
    :CommonSolve,
    :DSP,
    :DelimitedFiles,
    :Distributions,
    :Downloads,
    :FiniteDiff,
    :Glob,
    :IJulia,
    :Interpolations,
    :Latexify,
    :LinearAlgebra,
    :Metaheuristics,
    :MonteCarloMeasurements,
    :NLopt,
    :NaNStatistics,
    :OrdinaryDiffEqRosenbrock,
    :OrdinaryDiffEqTsit5,
    :Pkg,
    :Plots,
    :SciMLBase,
    :SciMLStructures,
    :Statistics,
    :StructArrays,
    :TOML,
    :Tables,
]

const RELEASE_PACKAGES = [
    :GlobalAnalysis,
    :Catalyst,
    :IJulia,
    :Latexify,
    :Metaheuristics,
    :MonteCarloMeasurements,
    :Plots,
]

"""
    packagesForKind(kind)

Returns the package list to include in the sysimage `kind`.
"""
function packagesForKind(kind::AbstractString)
    kind = gaNormalizeKind(kind)
    if kind == "dev"
        # Keep GlobalAnalysis out of the dev image so local source edits load
        # normally while the heavy dependencies still come from the sysimage.
        return DEV_PACKAGES
    else
        # Release/demo images are snapshots: maximum startup speed, but they
        # must be rebuilt when GlobalAnalysis source changes.
        return RELEASE_PACKAGES
    end
end

"""
    workloadForKind(kind)

Returns the precompile workload script path for the sysimage `kind`.
"""
function workloadForKind(kind::AbstractString)
    kind = gaNormalizeKind(kind)
    workload_file = kind == "dev" ? "workloadDev.jl" : "workloadRelease.jl"
    return joinpath(GA_SYSIMAGE_SCRIPT_DIR, workload_file)
end

"""
    sysimageBuildArgs()

Returns additional Julia arguments used while compiling the sysimage object.
"""
function sysimageBuildArgs()
    opt_level = strip(get(ENV, "GA_SYSIMAGE_OPT_LEVEL", ""))
    isempty(opt_level) && return ``

    opt_level in ("0", "1", "2", "3") ||
        throw(ArgumentError("GA_SYSIMAGE_OPT_LEVEL must be 0, 1, 2, or 3"))

    return Cmd(["-O$(opt_level)"])
end

"""
    hostCpuFeaturesBuildTarget(kind, cpu_target)

Returns the compile-time CPU target preference to give HostCPUFeatures while
building a sysimage, or `nothing` when no override is needed.
"""
function hostCpuFeaturesBuildTarget(kind::AbstractString, cpu_target::AbstractString)
    kind = gaNormalizeKind(kind)
    if kind != "release" || occursin("native", cpu_target)
        return nothing
    end

    return first(split(cpu_target, ';'))
end

"""
    withHostCpuFeaturesPreferences(kind, cpu_target) do
        ...
    end

Temporarily writes HostCPUFeatures compile-time preferences for portable release
sysimages. PackageCompiler runs dependency precompilation and tracing in helper
Julia processes that otherwise see the default runtime CPU target (`native`),
which can freeze host-specific HostCPUFeatures methods into the sysimage even
when the final object is emitted with a generic CPU target.
"""
function withHostCpuFeaturesPreferences(f::Function, kind::AbstractString, cpu_target::AbstractString)
    preference_target = hostCpuFeaturesBuildTarget(kind, cpu_target)
    preference_target === nothing && return f()

    project_path = joinpath(GA_REPO_ROOT, "Project.toml")
    preferences_path = joinpath(GA_REPO_ROOT, "LocalPreferences.toml")
    previous_project = read(project_path, String)
    had_preferences = isfile(preferences_path)
    previous_preferences = had_preferences ? read(preferences_path, String) : nothing

    project = TOML.parse(previous_project)
    extras = get!(project, "extras", Dict{String,Any}())
    extras[HOST_CPU_FEATURES_PREFS] = HOST_CPU_FEATURES_UUID
    open(project_path, "w") do io
        TOML.print(io, project; sorted=true)
    end

    preferences = had_preferences ? TOML.parse(previous_preferences) : Dict{String,Any}()
    host_preferences = get!(preferences, HOST_CPU_FEATURES_PREFS, Dict{String,Any}())
    host_preferences["cpu_target"] = preference_target
    host_preferences["freeze_cpu_target"] = true

    println("HostCPUFeatures CPU target preference: ", preference_target)
    open(preferences_path, "w") do io
        TOML.print(io, preferences; sorted=true)
    end

    try
        return f()
    finally
        write(project_path, previous_project)
        if had_preferences
            write(preferences_path, previous_preferences)
        else
            rm(preferences_path; force=true)
        end
    end
end

kind = gaKindFromArgs()
mkpath(GA_BUILD_DIR)

default_sysimage = gaDefaultSysimagePath(kind)
sysimage_path = get(ENV, "GA_SYSIMAGE_PATH", default_sysimage)
cpu_target = get(ENV, "JULIA_CPU_TARGET", gaDefaultCpuTarget(kind))
build_args = sysimageBuildArgs()
workload = workloadForKind(kind)
packages = packagesForKind(kind)

println("Repository: ", GA_REPO_ROOT)
println("Kind:       ", kind)
println("Workload:   ", workload)
println("Sysimage:   ", sysimage_path)
println("CPU target: ", cpu_target)
println("Build args: ", isempty(build_args.exec) ? "(default)" : join(build_args.exec, " "))
println("Packages:   ", join(string.(packages), ", "))

withHostCpuFeaturesPreferences(kind, cpu_target) do
    create_sysimage(
        packages;
        project=GA_REPO_ROOT,
        sysimage_path=sysimage_path,
        precompile_execution_file=workload,
        cpu_target=cpu_target,
        sysimage_build_args=build_args,
        incremental=true,
    )
end

println("Wrote sysimage: ", sysimage_path)
println("Run `julia --project=scripts/sysimage scripts/sysimage/stage.jl $(kind)` to install it in the default staged location.")

if gaNeedsWindowsStaging(abspath(sysimage_path))
    println("Windows note: this sysimage path contains spaces or non-ASCII characters.")
    println("Staging avoids fragile --sysimage paths in notebook kernels.")
end
