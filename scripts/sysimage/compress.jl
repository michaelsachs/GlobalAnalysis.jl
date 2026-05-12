include("common.jl")

kind = gaKindFromArgs("release")
sysimage_path = abspath(get(ENV, "GA_SYSIMAGE_PATH", gaDefaultSysimagePath(kind)))
compressed_path = abspath(get(ENV, "GA_SYSIMAGE_ARCHIVE_PATH", gaDefaultArchivePath(kind)))

isfile(sysimage_path) || error("Sysimage not found: $sysimage_path")
mkpath(dirname(compressed_path))
rm(compressed_path; force=true)

tar_cmd = `tar -czf $compressed_path -C $(dirname(sysimage_path)) $(basename(sysimage_path))`
if Sys.iswindows() && success(pipeline(`tar --force-local --help`; stdout=devnull, stderr=devnull))
    run(`tar --force-local -czf $compressed_path -C $(dirname(sysimage_path)) $(basename(sysimage_path))`)
else
    run(tar_cmd)
end

println("Wrote compressed sysimg: ", compressed_path)
println("Compressed size:        ", filesize(compressed_path))
