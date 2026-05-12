# staging a sysimg installs the release sysimg from build/sysimages 
# into the user sysimg folder so Julia/Jupyter can use it reliably 

include("common.jl")

kind = gaKindFromArgs()
source = get(ENV, "GA_SYSIMAGE_PATH", gaDefaultSysimagePath(kind))
source = abspath(source)

isfile(source) || error("Sysimage not found: $source")

stage_dir = get(ENV, "GA_SYSIMAGE_STAGE_DIR", gaDefaultStageDir())
mkpath(stage_dir)

staged = abspath(joinpath(stage_dir, basename(source)))

if staged != source
    rm(staged; force=true)

    try
        hardlink(source, staged)
        println(stderr, "Staged sysimage hardlink: $staged")
    catch err
        println(stderr, "Hardlink failed; copying sysimage instead.")
        showerror(stderr, err)
        println(stderr)
        cp(source, staged; force=true)
        println(stderr, "Staged sysimage copy: $staged")
    end
end

println(staged)
