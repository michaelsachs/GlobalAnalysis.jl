using Test
using GlobalAnalysis

@testset "Notebook Workflow Setup" begin
    autoThreads = GlobalAnalysis.getThreads("auto")

    @test autoThreads isa Int
    @test autoThreads > 0
    @test autoThreads == Sys.CPU_THREADS
    @test GlobalAnalysis.getThreads(:auto) == autoThreads
    @test GlobalAnalysis.getThreads(4) == 4
    @test GlobalAnalysis.getThreads("8") == 8

    @test_throws ArgumentError GlobalAnalysis.getThreads(0)
    @test_throws ArgumentError GlobalAnalysis.getThreads("many")

    notebookDir = mktempdir()
    sourceDir = mktempdir()
    write(joinpath(sourceDir, "example.ipynb"), "starter")
    write(joinpath(sourceDir, "notes.md"), "notes")

    copiedNotebooks = GlobalAnalysis.copyNotebooks(notebookDir; sourcePath=sourceDir)

    @test sort(basename.(copiedNotebooks)) == ["example.ipynb", "notes.md"]
    @test read(joinpath(notebookDir, "example.ipynb"), String) == "starter"

    write(joinpath(sourceDir, "example.ipynb"), "updated")
    copiedNotebooks = GlobalAnalysis.copyNotebooks(notebookDir; sourcePath=sourceDir)

    @test copiedNotebooks == String[]
    @test read(joinpath(notebookDir, "example.ipynb"), String) == "starter"

    @test GlobalAnalysis.defaultNotebookPath() == normpath(joinpath(homedir(), "GlobalAnalysis.jl"))
    @test GlobalAnalysis.resolveNotebookPath(notebookDir) == normpath(abspath(notebookDir))
    @test GlobalAnalysis.kernelName == "JuliaGA"
    @test GlobalAnalysis.kernelArguments() == ["--project=$(GlobalAnalysis.packageProjectPath())"]

    @test first(methods(GlobalAnalysis.ensureJupyterInstalled)).module == GlobalAnalysis
    @test first(methods(setup)).module == GlobalAnalysis
    @test first(methods(launch)).module == GlobalAnalysis
    @test_throws ArgumentError launch(notebookPath=joinpath(notebookDir, "missing"))
end
