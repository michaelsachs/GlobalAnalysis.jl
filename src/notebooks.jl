using Downloads
using IJulia
using Pkg
using TOML

const kernelName = "JuliaGA"
const jupyterCondaPackage = "notebook<7"
const sysimageReleaseBaseUrl = "https://github.com/michaelsachs/GlobalAnalysis.jl/releases/download"

"""
    packageNotebookPath()

Returns the notebook folder distributed with GlobalAnalysis.jl.
"""
packageNotebookPath() = normpath(joinpath(packageProjectPath(), "notebooks"))


"""
    defaultNotebookPath()

Returns the default user notebook folder.
"""
defaultNotebookPath() = normpath(joinpath(homedir(), "GlobalAnalysis.jl"))


"""
    packageProjectPath()

Returns the package project folder used when registering the notebook kernel.
"""
function packageProjectPath()
    packageFile = Base.find_package(string(nameof(@__MODULE__)))
    packageFile !== nothing && return normpath(dirname(dirname(packageFile)))

    return normpath(joinpath(@__DIR__, ".."))
end


"""
    sysimageLibraryExtension()

Returns the platform-specific shared library extension used by Julia
sysimages.
"""
function sysimageLibraryExtension()
    if Sys.iswindows()
        return ".dll"
    elseif Sys.isapple()
        return ".dylib"
    else
        return ".so"
    end
end


"""
    sysimagePlatformTag()

Returns the Julia version and platform tag used in GlobalAnalysis sysimage
filenames.
"""
function sysimagePlatformTag()
    os = Sys.iswindows() ? "windows" : Sys.isapple() ? "macos" : Sys.islinux() ? "linux" : Sys.KERNEL
    return "julia$(VERSION.major).$(VERSION.minor)-$(os)-$(Sys.ARCH)"
end


"""
    defaultSysimageStageDir()

Returns the directory where locally staged GlobalAnalysis sysimages are
expected.
"""
function defaultSysimageStageDir()
    if Sys.iswindows()
        return joinpath(get(ENV, "LOCALAPPDATA", tempdir()), "GlobalAnalysis.jl", "sysimages")
    elseif Sys.isapple()
        return joinpath(homedir(), "Library", "Application Support", "GlobalAnalysis.jl", "sysimages")
    else
        xdgDataHome = get(ENV, "XDG_DATA_HOME", joinpath(homedir(), ".local", "share"))
        return joinpath(xdgDataHome, "GlobalAnalysis.jl", "sysimages")
    end
end


"""
    defaultSysimagePath(kind=:dev)

Returns the default staged sysimage path for `kind`, which may be `:dev` or
`:release`.
"""
function defaultSysimagePath(kind::Union{AbstractString,Symbol}=:dev)
    kindString = lowercase(string(kind))
    kindString in ("dev", "release") ||
        throw(ArgumentError("sysimage kind must be :dev or :release"))

    filename = "GlobalAnalysis-$(kindString)-$(sysimagePlatformTag())$(sysimageLibraryExtension())"
    return normpath(joinpath(defaultSysimageStageDir(), filename))
end


"""
    defaultSysimageArchivePath(kind=:release)

Returns the default compressed archive path for the staged sysimage `kind`.
"""
function defaultSysimageArchivePath(kind::Union{AbstractString,Symbol}=:release)
    return defaultSysimagePath(kind) * ".tar.gz"
end


"""
    packageVersion()

Returns the package version from the GlobalAnalysis project file.
"""
function packageVersion()
    project = TOML.parsefile(joinpath(packageProjectPath(), "Project.toml"))
    return string(project["version"])
end


"""
    defaultSysimageReleaseTag()

Returns the release tag used when downloading prebuilt sysimage artifacts.
"""
function defaultSysimageReleaseTag()
    return get(ENV, "GA_SYSIMAGE_RELEASE_TAG", "v$(packageVersion())")
end


"""
    sysimageReleaseDownloadUrl(kind=:release; tag=defaultSysimageReleaseTag())

Returns the GitHub release asset URL for a compressed sysimage archive.
"""
function sysimageReleaseDownloadUrl(kind::Union{AbstractString,Symbol}=:release; tag=defaultSysimageReleaseTag())
    archiveName = basename(defaultSysimageArchivePath(kind))
    return "$(sysimageReleaseBaseUrl)/$(tag)/$(archiveName)"
end


"""
    extractSysimageArchive!(archivePath, sysimagePath)

Extracts a compressed sysimage archive into the staged sysimage directory.
"""
function extractSysimageArchive!(archivePath, sysimagePath)
    tar = Sys.which("tar")
    tar === nothing && error("Could not find `tar`; cannot extract sysimage archive")

    mkpath(dirname(sysimagePath))
    rm(sysimagePath; force=true)
    run(`$tar -xzf $archivePath -C $(dirname(sysimagePath))`)
    isfile(sysimagePath) || error("Archive did not contain expected sysimage: $sysimagePath")

    return sysimagePath
end


"""
    downloadSysimage(kind=:release; force=false, tag=defaultSysimageReleaseTag())

Downloads and extracts the compressed release sysimage for the current Julia
version and platform.
"""
function downloadSysimage(kind::Union{AbstractString,Symbol}=:release; force=false, tag=defaultSysimageReleaseTag())
    sysimagePath = defaultSysimagePath(kind)
    if isfile(sysimagePath) && !force
        return sysimagePath
    end

    archivePath = defaultSysimageArchivePath(kind)
    downloadPath = archivePath * ".download"
    url = sysimageReleaseDownloadUrl(kind; tag=tag)

    mkpath(dirname(archivePath))
    rm(downloadPath; force=true)

    try
        println("Downloading GlobalAnalysis sysimage from $(url)")
        Downloads.download(url, downloadPath)
        mv(downloadPath, archivePath; force=true)
        extractSysimageArchive!(archivePath, sysimagePath)
    finally
        rm(downloadPath; force=true)
    end

    return sysimagePath
end


"""
    kernelSysimageLabel(sysimage)

Returns the short label used in kernel names for the selected `sysimage`
mode.
"""
function kernelSysimageLabel(sysimage)
    if sysimage === nothing || sysimage === false
        return nothing
    elseif sysimage === true
        return "dev sysimg"
    elseif sysimage isa Symbol && sysimage === :dev
        return "dev sysimg"
    elseif sysimage isa Symbol && sysimage === :release
        return "sysimg"
    elseif sysimage isa AbstractString
        stripped = strip(sysimage)
        lowered = lowercase(stripped)
        if lowered == "dev"
            return "dev sysimg"
        elseif lowered == "release"
            return "sysimg"
        else
            return "sysimg"
        end
    else
        throw(ArgumentError("sysimage must be nothing, false, true, :dev, :release, or a path string"))
    end
end


"""
    kernelSpecName(; sysimage=nothing)

Returns the Jupyter kernelspec name used by the GlobalAnalysis notebook
kernel.
"""
function kernelSpecName(; sysimage=nothing)
    label = kernelSysimageLabel(sysimage)
    suffix = label === nothing ? "" : "-$(replace(label, " " => "-"))"
    return "$(lowercase(kernelName))-$(VERSION.major).$(VERSION.minor)$(suffix)"
end


"""
    kernelDisplayName(; sysimage=nothing)

Returns the user-facing Jupyter kernel name shown in the browser.
"""
function kernelDisplayName(; sysimage=nothing)
    label = kernelSysimageLabel(sysimage)
    if label === nothing
        return "$(kernelName) $(VERSION.major).$(VERSION.minor)"
    end

    return "$(kernelName) $(VERSION.major).$(VERSION.minor) ($(label))"
end


"""
    preparePackageProject()

Instantiates the package project used by the notebook kernel.
"""
function preparePackageProject()
    activeProject = Base.active_project()
    packageProject = packageProjectPath()

    try
        # setup may be called from a temporary outer project; prepare the
        # project that the installed Jupyter kernel will actually launch
        Pkg.activate(packageProject)
        Pkg.instantiate()
    finally
        if activeProject !== nothing
            Pkg.activate(activeProject)
        end
    end
end


"""
    getThreads(threads)

Returns the number of Julia threads to use. If `threads` is `"auto"` or
`:auto`, all available CPU threads are used.
"""
function getThreads(threads)

    errString = "threads must be \"auto\" or a positive integer"

    # threads entered as number
    if threads isa Integer
        if threads ≤ 0
            throw(ArgumentError(errString))
        end
        return threads
    end

    # threads entered as symbol
    if threads isa Symbol
        threads === :auto || throw(ArgumentError(errString))
        return Sys.CPU_THREADS
    end

    # threads entered as string
    if threads isa AbstractString
        threadsString = strip(threads)
        lowercase(threadsString) == "auto" && return Sys.CPU_THREADS

        threadsInt = tryparse(Int, threadsString)
        if threadsInt !== nothing && threadsInt > 0
            return threadsInt
        end
    end

    throw(ArgumentError(errString))
end


"""
    resolveNotebookPath(notebookPath)

Expands and normalizes `notebookPath` to an absolute path.
"""
function resolveNotebookPath(notebookPath)
    notebookPath isa AbstractString ||
        throw(ArgumentError("notebookPath must be a path string"))

    # path in user directory
    path = normpath(abspath(expanduser(notebookPath)))

    return path
end


"""
    makeFileWritable!(path)

Ensures a copied notebook can be updated after package installs that mark
source files as read-only.
"""
function makeFileWritable!(path)
    isfile(path) || return path

    if filemode(path) & 0o200 == 0
        chmod(path, filemode(path) | 0o200)
    end

    return path
end


"""
    copyNotebooks(notebookPath; sourcePath=packageNotebookPath())

Copies notebooks from `sourcePath` into `notebookPath`. Existing files are
left unchanged.
"""
function copyNotebooks(notebookPath; sourcePath=packageNotebookPath())
    notebookDir = resolveNotebookPath(notebookPath)
    sourceDir = resolveNotebookPath(sourcePath)
    isdir(notebookDir) || mkpath(notebookDir)

    copiedNotebooks = String[]
    for notebook in readdir(sourceDir)
        sourceFile = joinpath(sourceDir, notebook)
        targetFile = joinpath(notebookDir, notebook)

        # do not overwrite notebooks the user may have edited
        if !ispath(targetFile)
            cp(sourceFile, targetFile)
            makeFileWritable!(targetFile)
            push!(copiedNotebooks, targetFile)
        end
    end

    println("✓ Copied notebooks to $(notebookDir)")

    return copiedNotebooks
end


"""
    setNotebookKernel!(notebookPath; sysimage=nothing)

Sets the Jupyter kernel metadata of `notebookPath` to the dedicated
GlobalAnalysis kernel.
"""
function setNotebookKernel!(notebookPath; sysimage=nothing)
    notebookFile = resolveNotebookPath(notebookPath)
    notebook = read(notebookFile, String)

    kernelSpec = """
"kernelspec": {
   "display_name": "$(kernelDisplayName(; sysimage=sysimage))",
   "language": "julia",
   "name": "$(kernelSpecName(; sysimage=sysimage))"
  }"""

    # only update the notebook metadata; code cells and outputs stay untouched
    updatedNotebook = replace(
        notebook,
        r"\"kernelspec\"\s*:\s*\{.*?\n\s*\}"s => kernelSpec;
        count=1,
    )

    if updatedNotebook != notebook
        makeFileWritable!(notebookFile)
        write(notebookFile, updatedNotebook)
        return true
    end

    return false
end


"""
    setNotebookKernels!(notebookPath; sysimage=nothing)

Sets all copied notebook files in `notebookPath` to use the dedicated
GlobalAnalysis kernel.
"""
function setNotebookKernels!(notebookPath; sysimage=nothing)
    notebookDir = resolveNotebookPath(notebookPath)
    updatedNotebooks = String[]

    for notebook in readdir(notebookDir)
        notebookFile = joinpath(notebookDir, notebook)

        # notebook metadata decides which kernel Jupyter connects by default
        if isfile(notebookFile) && endswith(lowercase(notebook), ".ipynb")
            if setNotebookKernel!(notebookFile; sysimage=sysimage)
                push!(updatedNotebooks, notebookFile)
            end
        end
    end

    println("✓ Set notebooks to use kernel \"$(kernelDisplayName(; sysimage=sysimage))\"")

    return updatedNotebooks
end


"""
    resolveKernelSysimage(sysimage)

Resolves a sysimage mode or path to the concrete sysimage path used by the
notebook kernel.
"""
function resolveKernelSysimage(sysimage)
    if sysimage === nothing || sysimage === false
        return nothing
    elseif sysimage === true
        return defaultSysimagePath(:dev)
    elseif sysimage isa Symbol && sysimage in (:dev, :release)
        return defaultSysimagePath(sysimage)
    elseif sysimage isa AbstractString
        stripped = strip(sysimage)
        if lowercase(stripped) in ("dev", "release")
            return defaultSysimagePath(stripped)
        else
            return normpath(abspath(expanduser(stripped)))
        end
    else
        throw(ArgumentError("sysimage must be nothing, false, true, :dev, :release, or a path string"))
    end
end


"""
    kernelArguments(; sysimage=nothing)

Returns the Julia command line arguments used by the GlobalAnalysis notebook
kernel.
"""
function kernelArguments(; sysimage=nothing)
    args = ["--project=$(packageProjectPath())"]
    sysimagePath = resolveKernelSysimage(sysimage)

    if sysimagePath !== nothing
        isfile(sysimagePath) ||
            throw(ArgumentError("sysimage does not exist: $sysimagePath"))
        push!(args, "--sysimage=$(sysimagePath)")
    end

    return args
end


"""
    resolveSetupSysimage(sysimage; downloadReleaseSysimage=true)

Resolves the sysimage preference used by `setup` and downloads the release
sysimage when requested.
"""
function resolveSetupSysimage(sysimage; downloadReleaseSysimage=true)
    if sysimage === :auto || (sysimage isa AbstractString && lowercase(strip(sysimage)) == "auto")
        releasePath = defaultSysimagePath(:release)
        if !isfile(releasePath) && downloadReleaseSysimage
            try
                downloadSysimage(:release)
            catch err
                @warn "Could not install release sysimage; falling back to the regular JuliaGA kernel" exception = (err, catch_backtrace())
            end
        end

        return isfile(releasePath) ? :release : nothing
    elseif sysimage isa Symbol && sysimage === :release
        releasePath = defaultSysimagePath(:release)
        if !isfile(releasePath) && downloadReleaseSysimage
            downloadSysimage(:release)
        end
    elseif sysimage isa AbstractString && lowercase(strip(sysimage)) == "release"
        releasePath = defaultSysimagePath(:release)
        if !isfile(releasePath) && downloadReleaseSysimage
            downloadSysimage(:release)
        end
    end

    return sysimage
end


"""
    installNotebookKernel(kernelEnv; sysimage=nothing)

Installs or updates a GlobalAnalysis Jupyter kernel using `kernelEnv` and the
selected sysimage mode.
"""
function installNotebookKernel(kernelEnv; sysimage=nothing)
    IJulia.installkernel(
        kernelName,
        kernelArguments(; sysimage=sysimage)...;
        specname=kernelSpecName(; sysimage=sysimage),
        displayname=kernelDisplayName(; sysimage=sysimage),
        env=kernelEnv,
    )

    println("✓ Added kernel \"$(kernelDisplayName(; sysimage=sysimage))\" to Jupyter ")

    return nothing
end


"""
    kernelEnvironment(threads)

Returns the environment variables used by the GlobalAnalysis notebook
kernel. `threads` accepts the same values as [`getThreads`](@ref), but is
normalized to a concrete integer before being written to the kernelspec.
"""
function kernelEnvironment(threads)
    threadsInt = getThreads(threads)

    return Dict(
        "JULIA_NUM_THREADS" => string(threadsInt),
        "JULIA_DEPOT_PATH" => join(DEPOT_PATH, Sys.iswindows() ? ";" : ":"),
    )
end


"""
    ensureJupyterInstalled()

Checks whether Jupyter is available to IJulia. If IJulia is using its
managed Conda environment and Jupyter is missing, it is installed.
"""
function ensureJupyterInstalled()
    notebookCmd = IJulia.find_jupyter_subcommand("notebook")
    jupyter = first(notebookCmd)

    scriptDir = IJulia.get_Conda() do Conda
        Conda.SCRIPTDIR
    end

    if dirname(jupyter) == abspath(scriptDir) && !Sys.isexecutable(IJulia.exe(jupyter, "-notebook"))
        # install once during setup so launch can open without prompting
        IJulia.get_Conda() do Conda
            Conda.add(jupyterCondaPackage)
        end
        
        println("✓ Installed Jupyter")
    else
        println("✓ Jupyter is already installed")
    end

end


"""
    setup(; threads="auto", notebookPath=defaultNotebookPath(), sysimage=:auto, downloadReleaseSysimage=true)

Registers the `JuliaGA` Jupyter kernel for browser notebooks.

`threads` controls the number of Julia threads used by the kernel and may be
`"auto"` or a positive integer. `"auto"` uses all available CPU threads.
`notebookPath` is created if it does not already exist.
`sysimage` defaults to `:auto`, which uses the release sysimage when it can be
installed and otherwise falls back to the regular kernel. It may also be
`:dev`, `:release`, `false`, `true` for the dev sysimage, or a path.
"""
function setup(; threads="auto", notebookPath=defaultNotebookPath(), sysimage=:auto, downloadReleaseSysimage=true)
    # determine thread count before writing the kernelspec
    kernelEnv = kernelEnvironment(threads)

    # prepare the package environment before the browser kernel starts
    preparePackageProject()

    # install jupyter if needed
    ensureJupyterInstalled()

    # keep a regular kernel available even when the sysimage kernel is default
    installNotebookKernel(kernelEnv; sysimage=nothing)
    defaultSysimage = resolveSetupSysimage(sysimage; downloadReleaseSysimage=downloadReleaseSysimage)
    if defaultSysimage !== nothing && defaultSysimage !== false
        installNotebookKernel(kernelEnv; sysimage=defaultSysimage)
    end

    notebookDir = resolveNotebookPath(notebookPath)
    copyNotebooks(notebookDir)

    dataDir = joinpath(notebookDir, "data")
    isdir(dataDir) || mkpath(dataDir)
    sourceDataDir = normpath(joinpath(packageProjectPath(), "data"))
    for dataFile in readdir(sourceDataDir)
        sourceFile = joinpath(sourceDataDir, dataFile)
        targetFile = joinpath(dataDir, dataFile)

        # keep any user-edited or user-added data files intact
        if isfile(sourceFile) && !ispath(targetFile)
            cp(sourceFile, targetFile)
        end
    end
    println("✓ Copied example data to $(dataDir)")

    setNotebookKernels!(notebookDir; sysimage=defaultSysimage)

    println("✓ Setup complete!")

    return nothing
end

"""
    launch(; notebookPath=defaultNotebookPath(), sysimage=:auto)

Opens the browser Jupyter notebook interface in `notebookPath`.
`sysimage=:auto` uses the staged release sysimage kernel if it is available.
"""
function launch(; notebookPath=defaultNotebookPath(), sysimage=:auto)
    notebookDir = resolveNotebookPath(notebookPath)
    isdir(notebookDir) || throw(ArgumentError("notebookPath does not exist: $notebookDir"))

    defaultSysimage = resolveSetupSysimage(sysimage; downloadReleaseSysimage=false)
    setNotebookKernels!(notebookDir; sysimage=defaultSysimage)

    return IJulia.notebook(dir=notebookDir)
end
