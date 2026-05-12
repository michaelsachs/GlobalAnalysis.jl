const GA_SYSIMAGE_SCRIPT_DIR = @__DIR__
const GA_REPO_ROOT = normpath(joinpath(GA_SYSIMAGE_SCRIPT_DIR, "..", ".."))
const GA_BUILD_DIR = joinpath(GA_REPO_ROOT, "build", "sysimages")

const GA_SYSIMAGE_KINDS = ("dev", "release")

"""
    gaNormalizeKind(kind)

Normalizes and validates the sysimage `kind`.
"""
function gaNormalizeKind(kind::AbstractString)
    normalized = lowercase(strip(kind))
    normalized in GA_SYSIMAGE_KINDS ||
        throw(ArgumentError("sysimage kind must be one of $(join(GA_SYSIMAGE_KINDS, ", "))"))
    return normalized
end

"""
    gaKindFromArgs(default="dev")

Returns the sysimage kind from command line arguments or `GA_SYSIMAGE_KIND`.
"""
function gaKindFromArgs(default::AbstractString="dev")
    return gaNormalizeKind(!isempty(ARGS) ? ARGS[1] : get(ENV, "GA_SYSIMAGE_KIND", default))
end

"""
    gaLibraryExtension()

Returns the platform-specific shared library extension used by Julia
sysimages.
"""
function gaLibraryExtension()
    if Sys.iswindows()
        return ".dll"
    elseif Sys.isapple()
        return ".dylib"
    else
        return ".so"
    end
end

"""
    gaPlatformTag()

Returns the Julia version and platform tag used in sysimage filenames.
"""
function gaPlatformTag()
    os = Sys.iswindows() ? "windows" : Sys.isapple() ? "macos" : Sys.islinux() ? "linux" : Sys.KERNEL
    return "julia$(VERSION.major).$(VERSION.minor)-$(os)-$(Sys.ARCH)"
end

"""
    gaSysimageBasename(kind)

Returns the platform-specific filename for the sysimage `kind`.
"""
function gaSysimageBasename(kind::AbstractString)
    kind = gaNormalizeKind(kind)
    return "GlobalAnalysis-$(kind)-$(gaPlatformTag())$(gaLibraryExtension())"
end

"""
    gaDefaultSysimagePath(kind)

Returns the default build output path for the sysimage `kind`.
"""
function gaDefaultSysimagePath(kind::AbstractString)
    return joinpath(GA_BUILD_DIR, gaSysimageBasename(kind))
end

"""
    gaDefaultArchivePath(kind)

Returns the default compressed archive path for the sysimage `kind`.
"""
function gaDefaultArchivePath(kind::AbstractString)
    return gaDefaultSysimagePath(kind) * ".tar.gz"
end

"""
    gaDefaultCpuTarget(kind)

Returns the default CPU target for the sysimage `kind`.
"""
function gaDefaultCpuTarget(kind::AbstractString)
    kind = gaNormalizeKind(kind)
    # Local dev images can be machine-specific. Release artifacts should run
    # across the runner architecture rather than only the exact CI CPU model.
    return kind == "dev" ? "native" : "generic"
end

"""
    gaPathIsAscii(path)

Returns whether `path` contains only ASCII characters.
"""
function gaPathIsAscii(path::AbstractString)
    return all(c -> Int(c) <= 0x7f, path)
end

"""
    gaNeedsWindowsStaging(path)

Returns whether `path` should be staged before use as a Windows sysimage path.
"""
function gaNeedsWindowsStaging(path::AbstractString)
    return Sys.iswindows() && (!gaPathIsAscii(path) || occursin(" ", path))
end

"""
    gaDefaultStageDir()

Returns the default user directory for staged GlobalAnalysis sysimages.
"""
function gaDefaultStageDir()
    if Sys.iswindows()
        return joinpath(get(ENV, "LOCALAPPDATA", tempdir()), "GlobalAnalysis.jl", "sysimages")
    elseif Sys.isapple()
        return joinpath(homedir(), "Library", "Application Support", "GlobalAnalysis.jl", "sysimages")
    else
        xdg_data_home = get(ENV, "XDG_DATA_HOME", joinpath(homedir(), ".local", "share"))
        return joinpath(xdg_data_home, "GlobalAnalysis.jl", "sysimages")
    end
end

"""
    gaStagedSysimagePath(kind)

Returns the default staged sysimage path for `kind`.
"""
function gaStagedSysimagePath(kind::AbstractString)
    return joinpath(gaDefaultStageDir(), gaSysimageBasename(kind))
end
