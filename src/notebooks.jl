using IJulia

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
            push!(copiedNotebooks, targetFile)
        end
    end

    println("✓ Copied notebooks to $(notebookDir)")
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

    # install jupyter if needed
    ensureJupyterInstalled()

    # env for kernel installation
    kernelEnv = Dict("JULIA_NUM_THREADS" => string(threadsInt))

    # install/update dedicated kernel for browser notebooks
    IJulia.installkernel(
        kernelName,
        kernelArguments()...;
        env=kernelEnv,
    )

    println("✓ Added kernel \"$(kernelName)\" to Jupyter ")

    notebookDir = resolveNotebookPath(notebookPath)
    copyNotebooks(notebookDir)

end

"""
    launch(; notebookPath=defaultNotebookPath())

Opens the browser Jupyter notebook interface in `notebookPath`.
"""
function launch(; notebookPath=defaultNotebookPath())
    notebookDir = resolveNotebookPath(notebookPath)
    isdir(notebookDir) || throw(ArgumentError("notebookPath does not exist: $notebookDir"))

    return IJulia.notebook(dir=notebookDir)
end
