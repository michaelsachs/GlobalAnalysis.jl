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
    gaPortableReleaseBuild(kind, cpu_target)

Returns whether this sysimage build is a release artifact with a non-native CPU
target.
"""
function gaPortableReleaseBuild(kind::AbstractString, cpu_target::AbstractString)
    kind = gaNormalizeKind(kind)
    return kind == "release" && !occursin("native", cpu_target)
end

"""
    ensureBuildProcessCpuTarget(kind, cpu_target)

Checks that release sysimage builds were started with the same CPU target that
will be passed to PackageCompiler. This catches the old local invocation style,
where only the final sysimage object used the portable CPU target.
"""
function ensureBuildProcessCpuTarget(kind::AbstractString, cpu_target::AbstractString)
    gaPortableReleaseBuild(kind, cpu_target) || return nothing

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
    gaPortableReleaseBuild(kind, cpu_target) || return nothing

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

"""
    withHostCpuFeaturesPreferences(kind, cpu_target) do
        ...
    end

Temporarily freezes HostCPUFeatures to the release CPU target. Packages such as
VectorizationBase generate method bodies from HostCPUFeatures constants; without
this, a cached host-specific feature set can emit instructions the portable
sysimage target cannot lower.
"""
function withHostCpuFeaturesPreferences(f::Function, kind::AbstractString, cpu_target::AbstractString)
    gaPortableReleaseBuild(kind, cpu_target) || return f()

    project_path = joinpath(GA_REPO_ROOT, "Project.toml")
    preferences_path = joinpath(GA_REPO_ROOT, "LocalPreferences.toml")
    previous_project = read(project_path, String)
    had_preferences = isfile(preferences_path)
    previous_preferences = had_preferences ? read(preferences_path, String) : nothing
    preference_target = first(split(cpu_target, ';'))

    try
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
        open(preferences_path, "w") do io
            TOML.print(io, preferences; sorted=true)
        end

        println("HostCPUFeatures CPU target preference: ", preference_target)
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

function hostCpuFeaturesAarch64Source()
    previous_load_path = copy(LOAD_PATH)
    package_file = try
        pushfirst!(LOAD_PATH, GA_REPO_ROOT)
        Base.locate_package(Base.PkgId(Base.UUID(HOST_CPU_FEATURES_UUID), HOST_CPU_FEATURES_PREFS))
    finally
        empty!(LOAD_PATH)
        append!(LOAD_PATH, previous_load_path)
    end

    package_file === nothing && return nothing
    return joinpath(dirname(package_file), "cpu_info_aarch64.jl")
end

function patchHostCpuFeaturesAarch64Source!(source::AbstractString)
    previous_source = read(source, String)
    patched_source = replace(
        previous_source,
        """@noinline vscale() = ccall("llvm.vscale.i64", llvmcall, Int64, ())""" =>
            """@noinline vscale() = one(Int64)""",
        """@noinline vscale() = ccall("llvm.vscale.i32", llvmcall, Int32, ())""" =>
            """@noinline vscale() = one(Int32)""",
    )

    patched_source == previous_source &&
        error("HostCPUFeatures AArch64 vscale source did not match expected text: $source")

    write(source, patched_source)
    return nothing
end

"""
    withMacHostCpuFeaturesVscaleSourcePatch(kind, cpu_target) do
        ...
    end

Temporarily rewrites HostCPUFeatures.vscale on macOS/AArch64 portable release
builds. The method is unused for Apple Silicon without SVE, but its unconditional
`llvm.vscale` definition is compiled while PackageCompiler writes the sysimage
object and LLVM aborts for the `apple-m1` target.
"""
function withMacHostCpuFeaturesVscaleSourcePatch(
    f::Function,
    kind::AbstractString,
    cpu_target::AbstractString,
)
    gaPortableReleaseBuild(kind, cpu_target) || return f()
    (Sys.isapple() && Sys.ARCH === :aarch64) || return f()

    source = hostCpuFeaturesAarch64Source()
    source === nothing && error("Could not locate HostCPUFeatures AArch64 source")

    package_dir = dirname(dirname(source))
    patched_package_dir = joinpath(GA_BUILD_DIR, "HostCPUFeatures-vscale-patch")
    patched_source = joinpath(patched_package_dir, "src", "cpu_info_aarch64.jl")
    manifest_path = joinpath(GA_REPO_ROOT, "Manifest.toml")
    previous_manifest = read(manifest_path, String)

    try
        rm(patched_package_dir; recursive=true, force=true)
        cp(package_dir, patched_package_dir)
        chmod(patched_source, 0o644)
        patchHostCpuFeaturesAarch64Source!(patched_source)

        manifest = TOML.parse(previous_manifest)
        host_entries = manifest["deps"][HOST_CPU_FEATURES_PREFS]
        length(host_entries) == 1 ||
            error("Expected one HostCPUFeatures manifest entry, found $(length(host_entries))")
        host_entry = only(host_entries)
        delete!(host_entry, "git-tree-sha1")
        host_entry["path"] = patched_package_dir
        open(manifest_path, "w") do io
            TOML.print(io, manifest; sorted=true)
        end

        println("Using patched HostCPUFeatures source for macOS AArch64 sysimage: ", patched_package_dir)
        return f()
    finally
        write(manifest_path, previous_manifest)
        rm(patched_package_dir; recursive=true, force=true)
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

withHostCpuFeaturesPreferences(kind, cpu_target) do
    withMacHostCpuFeaturesVscaleSourcePatch(kind, cpu_target) do
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
end

println("Wrote sysimage: ", sysimage_path)
println("Run `julia --project=scripts/sysimage scripts/sysimage/stage.jl $(kind)` to install it in the default staged location.")

if gaNeedsWindowsStaging(abspath(sysimage_path))
    println("Windows note: this sysimage path contains spaces or non-ASCII characters.")
    println("Staging avoids fragile --sysimage paths in notebook kernels.")
end
