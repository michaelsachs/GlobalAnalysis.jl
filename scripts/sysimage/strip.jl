# stripping a sysimg removes nonessential binary metadata from the compiled sysimg file
# (debug symbols, symbol table entries not needed at runtime, linker/debug metadata, etc.).
# Currently saves ~120 MB before compression.

include("common.jl")

"""
    findStripTool()

Returns the available strip executable used to reduce sysimage size.
"""
function findStripTool()
    names = Sys.iswindows() ? ("strip.exe", "llvm-strip.exe") : ("strip", "llvm-strip")
    for name in names
        tool = Sys.which(name)
        tool === nothing || return tool
    end

    # PackageCompiler downloads a MinGW artifact on Windows. It is not always
    # on PATH, so look in the common artifact layout without depending on Glob.
    if Sys.iswindows()
        for depot in DEPOT_PATH
            artifacts = joinpath(depot, "artifacts")
            isdir(artifacts) || continue
            for hash in readdir(artifacts)
                candidate = joinpath(artifacts, hash, "extracted_files", "mingw64", "bin", "strip.exe")
                isfile(candidate) && return candidate
            end
        end
    end

    return nothing
end

kind = gaKindFromArgs("release")
sysimage_path = abspath(get(ENV, "GA_SYSIMAGE_PATH", gaDefaultSysimagePath(kind)))
isfile(sysimage_path) || error("Sysimage not found: $sysimage_path")

tool = findStripTool()
if tool === nothing
    println("No strip tool found; leaving sysimage unstripped: ", sysimage_path)
    exit()
end

before = filesize(sysimage_path)
if Sys.isapple()
    run(`$tool -x $sysimage_path`)
else
    run(`$tool --strip-unneeded $sysimage_path`)
end
after = filesize(sysimage_path)

println("Stripped sysimage: ", sysimage_path)
println("Size before: ", before)
println("Size after:  ", after)
