using IJulia
using Pkg

const kernelName = "JuliaGA"

"""
    packageNotebookPath()

Returns the notebook folder distributed with GlobalAnalysis.jl.
"""
packageNotebookPath() = normpath(joinpath(@__DIR__, "..", "notebooks"))


"""
    defaultNotebookPath()

Returns the default user notebook folder.
"""
defaultNotebookPath() = normpath(joinpath(homedir(), "GlobalAnalysis.jl"))


"""
    packageProjectPath()

Returns the package project folder used when registering the notebook kernel.
"""
packageProjectPath() = normpath(joinpath(@__DIR__, ".."))


"""
    kernelSpecName()

Returns the Jupyter kernelspec name used by the GlobalAnalysis notebook
kernel.
"""
kernelSpecName() = "$(lowercase(kernelName))-$(VERSION.major).$(VERSION.minor)"


"""
    kernelDisplayName()

Returns the user-facing Jupyter kernel name shown in the browser.
"""
kernelDisplayName() = "$(kernelName) $(VERSION.major).$(VERSION.minor)"


"""
    preparePackageProject()

Instantiates and precompiles the package project used by the notebook kernel.
"""
function preparePackageProject()
    activeProject = Base.active_project()
    packageProject = packageProjectPath()

    try
        # setup may be called from a temporary outer project; prepare the
        # project that the installed Jupyter kernel will actually launch
        Pkg.activate(packageProject)
        Pkg.instantiate()
        Pkg.precompile()
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
    setNotebookKernel!(notebookPath)

Sets the Jupyter kernel metadata of `notebookPath` to the dedicated
GlobalAnalysis kernel.
"""
function setNotebookKernel!(notebookPath)
    notebookFile = resolveNotebookPath(notebookPath)
    notebook = read(notebookFile, String)

    kernelSpec = """
"kernelspec": {
   "display_name": "$(kernelDisplayName())",
   "language": "julia",
   "name": "$(kernelSpecName())"
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
    setNotebookKernels!(notebookPath)

Sets all copied notebook files in `notebookPath` to use the dedicated
GlobalAnalysis kernel.
"""
function setNotebookKernels!(notebookPath)
    notebookDir = resolveNotebookPath(notebookPath)
    updatedNotebooks = String[]

    for notebook in readdir(notebookDir)
        notebookFile = joinpath(notebookDir, notebook)

        # notebook metadata decides which kernel Jupyter connects by default
        if isfile(notebookFile) && endswith(lowercase(notebook), ".ipynb")
            if setNotebookKernel!(notebookFile)
                push!(updatedNotebooks, notebookFile)
            end
        end
    end

    println("✓ Set notebooks to use kernel \"$(kernelDisplayName())\"")

    return updatedNotebooks
end


"""
    kernelArguments()

Returns the Julia command line arguments used by the GlobalAnalysis notebook
kernel.
"""
kernelArguments() = ["--project=$(packageProjectPath())"]


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
            Conda.add("jupyter")
        end
        
        println("✓ Installed Jupyter")
    else
        println("✓ Jupyter is already installed")
    end

end


"""
    setup(; threads="auto", notebookPath=defaultNotebookPath())

Registers the `JuliaGA` Jupyter kernel for browser notebooks.

`threads` controls the number of Julia threads used by the kernel and may be
`"auto"` or a positive integer. `"auto"` uses all available CPU threads.
`notebookPath` is created if it does not already exist.
"""
function setup(; threads="auto", notebookPath=defaultNotebookPath())
    # determine thread count before writing the kernelspec
    threadsInt = getThreads(threads)

    # prepare the package environment before the browser kernel starts
    preparePackageProject()

    # install jupyter if needed
    ensureJupyterInstalled()

    # env for kernel installation
    kernelEnv = Dict(
        "JULIA_NUM_THREADS" => string(threads),
        "JULIA_DEPOT_PATH" => join(DEPOT_PATH, Sys.iswindows() ? ";" : ":"),
    )

    # install/update dedicated kernel for browser notebooks
    IJulia.installkernel(
        kernelName,
        kernelArguments()...;
        env=kernelEnv,
    )

    println("✓ Added kernel \"$(kernelName)\" to Jupyter ")

    notebookDir = resolveNotebookPath(notebookPath)
    copyNotebooks(notebookDir)

    dataDir = joinpath(notebookDir, "data")
    isdir(dataDir) || mkpath(dataDir)
    sourceDataDir = normpath(joinpath(@__DIR__, "..", "data"))
    for dataFile in readdir(sourceDataDir)
        sourceFile = joinpath(sourceDataDir, dataFile)
        targetFile = joinpath(dataDir, dataFile)

        # keep any user-edited or user-added data files intact
        if isfile(sourceFile) && !ispath(targetFile)
            cp(sourceFile, targetFile)
        end
    end
    println("✓ Copied example data to $(dataDir)")

    setNotebookKernels!(notebookDir)

    println("✓ Setup complete!")

    return nothing
end

"""
    launch(; notebookPath=defaultNotebookPath())

Opens the browser Jupyter notebook interface in `notebookPath`.
"""
function launch(; notebookPath=defaultNotebookPath())
    notebookDir = resolveNotebookPath(notebookPath)
    isdir(notebookDir) || throw(ArgumentError("notebookPath does not exist: $notebookDir"))

    setNotebookKernels!(notebookDir)

    return IJulia.notebook(dir=notebookDir)
end
