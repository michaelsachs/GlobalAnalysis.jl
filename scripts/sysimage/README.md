# Sysimage Workflow

This directory supports two PackageCompiler sysimage modes.

`dev` sysimages include heavy dependencies but deliberately exclude
`GlobalAnalysis`, so local edits to package source load normally.

`release` sysimages include `GlobalAnalysis` and trace a notebook-like
workload for maximum startup speed. Rebuild release images after changing
package source.

## Local Dev Build

```powershell
julia --project=scripts/sysimage scripts/sysimage/build.jl dev
julia --project=scripts/sysimage scripts/sysimage/stage.jl dev
```

The staging step installs the image into the default user sysimage directory
and prints the path that notebook kernels pass to `--sysimage`.

For browser notebooks, install a kernel that uses the dev sysimage:

```julia
using GlobalAnalysis
setup(sysimage=:dev)
```

This installs a separate kernel such as `JuliaGA 1.12 (dev sysimg)`, so the
normal `JuliaGA 1.12` kernel remains available for comparison or fallback.
Plain `setup()` registers the normal kernel and then tries to install the
`JuliaGA 1.12 (sysimg)` kernel as the default notebook kernel.

## Release Build

```powershell
$env:JULIA_CPU_TARGET = "generic"
julia --cpu-target="$env:JULIA_CPU_TARGET" --project=scripts/sysimage scripts/sysimage/build.jl release
julia --project=scripts/sysimage scripts/sysimage/strip.jl release
julia --project=scripts/sysimage scripts/sysimage/compress.jl release
```

On Apple Silicon release builders, use `apple-m1` instead of `generic`. The
build pins PackageCompiler helper processes and HostCPUFeatures preferences to
the same CPU target. It also temporarily uses a patched local copy of
HostCPUFeatures while compiling the sysimage so LLVM does not lower the package's
unused AArch64 `llvm.vscale` method for Apple Silicon targets.

The compression step writes a `.tar.gz` file. GitHub Actions uploads that
compressed sysimg as a workflow artifact and, for published GitHub releases,
attaches it to the release so `setup()` can download it.
