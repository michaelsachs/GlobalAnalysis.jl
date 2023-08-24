using IJulia

function jupyterNotebook()
    dir = joinpath(@__DIR__, "..", "notebooks")
    notebook(dir=dir)
end

function jupyterLab()
    dir = joinpath(@__DIR__, "..", "notebooks")
    jupyterlab(dir=dir)
end