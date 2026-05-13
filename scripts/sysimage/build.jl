using PackageCompiler

include("common.jl")

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
    ensureBuildProcessCpuTarget(kind, cpu_target)

Checks that release sysimage builds were started with the same CPU target that
will be passed to PackageCompiler. This catches the old local invocation style,
where only the final sysimage object used the portable CPU target.
"""
function ensureBuildProcessCpuTarget(kind::AbstractString, cpu_target::AbstractString)
    kind = gaNormalizeKind(kind)
    (kind == "release" && !occursin("native", cpu_target)) || return nothing

    process_target = Base.unsafe_string(Base.JLOptions().cpu_target)
    process_target == cpu_target && return nothing

    throw(ArgumentError(
        "release sysimage builds must start Julia with --cpu-target=$(repr(cpu_target)); " *
        "current process CPU target is $(repr(process_target))"
    ))
end

"""
    alignPackageCompilerHelperCpuTarget!(kind, cpu_target)

Makes PackageCompiler's helper Julia processes use the same CPU target as the
final sysimage object. PackageCompiler passes `cpu_target` to the final object
compiler, but its package-precompile and workload-tracing helpers otherwise
start with the default target and can record host-only CPU feature methods.
"""
function alignPackageCompilerHelperCpuTarget!(kind::AbstractString, cpu_target::AbstractString)
    kind = gaNormalizeKind(kind)
    (kind == "release" && !occursin("native", cpu_target)) || return nothing

    @eval PackageCompiler begin
        _GA_HELPER_CPU_TARGET = $cpu_target

        function get_julia_cmd()
            julia_path = joinpath(Sys.BINDIR, Base.julia_exename())
            color = Base.have_color === nothing ? "auto" : Base.have_color ? "yes" : "no"
            cmd = if isdefined(Base, :Linking)
                `$julia_path --color=$color --startup-file=no --pkgimages=no`
            else
                `$julia_path --color=$color --startup-file=no`
            end
            return `$cmd --cpu-target=$_GA_HELPER_CPU_TARGET`
        end
    end

    return nothing
end

kind = gaKindFromArgs()
mkpath(GA_BUILD_DIR)

default_sysimage = gaDefaultSysimagePath(kind)
sysimage_path = get(ENV, "GA_SYSIMAGE_PATH", default_sysimage)
cpu_target = get(ENV, "JULIA_CPU_TARGET", gaDefaultCpuTarget(kind))
build_args = sysimageBuildArgs()
workload = workloadForKind(kind)
packages = packagesForKind(kind)

ensureBuildProcessCpuTarget(kind, cpu_target)
alignPackageCompilerHelperCpuTarget!(kind, cpu_target)

println("Repository: ", GA_REPO_ROOT)
println("Kind:       ", kind)
println("Workload:   ", workload)
println("Sysimage:   ", sysimage_path)
println("CPU target: ", cpu_target)
println("Process CPU target: ", Base.unsafe_string(Base.JLOptions().cpu_target))
println("Build args: ", isempty(build_args.exec) ? "(default)" : join(build_args.exec, " "))
println("Packages:   ", join(string.(packages), ", "))

create_sysimage(
    packages;
    project=GA_REPO_ROOT,
    sysimage_path=sysimage_path,
    precompile_execution_file=workload,
    cpu_target=cpu_target,
    sysimage_build_args=build_args,
    incremental=true,
)

println("Wrote sysimage: ", sysimage_path)
println("Run `julia --project=scripts/sysimage scripts/sysimage/stage.jl $(kind)` to install it in the default staged location.")

if gaNeedsWindowsStaging(abspath(sysimage_path))
    println("Windows note: this sysimage path contains spaces or non-ASCII characters.")
    println("Staging avoids fragile --sysimage paths in notebook kernels.")
end
