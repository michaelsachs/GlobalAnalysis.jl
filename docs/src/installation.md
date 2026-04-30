# Installation

This package is designed to take advantage of all available CPU cores to speed up fits and processing. Julia calls the parallel workers it uses for this `threads`. This package performs best if the number of Julia threads is set to the number of physical CPU cores.

The number of physical CPU cores is not straightforward to determine automatically, so it is best to provide it manually during setup. To figure out how many cores your CPU has, search for something like *"How many cores does my Intel i5-10400 have?"*. If you do not provide the number of cores, the setup will figure out a number automatically; if that number is different from the number of physical cores, you may get slightly slower fits, but the results will remain the same.

## Installing Julia

Download and install the current stable Julia release from <https://julialang.org/downloads/>. The recommended installer is Juliaup, which is the default installation method on the Julia download page.

After installation, open a terminal:

- Windows: open PowerShell or Windows Terminal.
- macOS: open Terminal.
- Linux: open your usual terminal.

Start Julia with the number of physical CPU cores you found. For example, for a CPU with 6 cores:

```bash
julia -t 6
```

If you do not know the number of physical CPU cores, use automatic threading instead:

```bash
julia -t auto
```

You should see a prompt that looks like this:

```julia
julia>
```

Starting Julia from the terminal is useful because the number of threads has to be set when Julia starts. `setup()` configures the future notebook kernel, but it cannot change the number of threads in the Julia session that is already running `setup()`.

If `julia -t 6` or `julia -t auto` is not found, close and reopen the terminal first. If it still does not work, Julia was not added to your system path during installation. In that case, you can still start Julia from the Julia application and continue with the steps below; setup may just take longer the first time.

## Installing GlobalAnalysis.jl

At the threaded Julia prompt, press `]` to enter Julia's package manager. The prompt changes from `julia>` to something ending in `pkg>`.

Then run:

```julia-repl
pkg> add https://github.com/michaelsachs/GlobalAnalysis.jl
```

This downloads GlobalAnalysis.jl and its Julia dependencies. The first installation can take a few minutes.

Press Backspace to return to the normal `julia>` prompt.

## Setting up notebooks

At the Julia prompt, first load the package:

```julia
using GlobalAnalysis
```

An initial one-time setup prepares the environment for the analysis notebooks. Before running the setup, decide whether you want to change the defaults for the browser notebook workflow:

- `threads`: number of Julia threads used by the notebook kernel. Choose the number of CPU cores if you know that number.
- `notebookPath`: folder where the example notebooks and example data are copied. If not provided, a `GlobalAnalysis.jl` folder will be created in your home directory and used for that purpose.

If you know the number of physical CPU cores, provide it to `setup()`. For example, for a CPU with 6 cores:

```julia
setup(threads=6)
```

To place the notebooks somewhere else, use:

```julia
setup(notebookPath=raw"C:\Users\YourName\Documents\GlobalAnalysis")
```

If the folder does not exist, it will be created.

You can also set both options at once:

```julia
setup(threads=6, notebookPath=raw"C:\Users\YourName\Documents\GlobalAnalysis")
```

Alternatively, if you want to run with the default options, simply use:

```julia
setup()
```

You can rerun `setup()` later if you want to change the thread count or notebook folder. Existing notebooks and data files are not overwritten, so your edits are safe when rerunning `setup()`.

## Launching notebooks

After setup, start the browser notebook interface. If you used a custom notebook folder, launch the same folder with:

```julia
launch(notebookPath=raw"C:\Users\YourName\Documents\GlobalAnalysis")
```

If not, just run:

```julia
launch()
```

Jupyter should open in your browser. The copied notebooks are configured to use the multithreaded `JuliaGA` kernel automatically. To get started, open `kineticModel.ipynb` and run the cells.

The first execution after starting a kernel will take longer because Julia compiles the required functions. Later executions in the same session are much faster.

## Updating

To update an installed version of GlobalAnalysis.jl, press `]` to enter the package manager and run:

```julia-repl
pkg> update GlobalAnalysis
```

Press Backspace to return to the normal `julia>` prompt.

Then restart the notebook kernel so the notebook loads the updated package code.

## Optional: Using VS Code for development

This step is not required for normal use. VS Code is useful if you want to edit the code of GlobalAnalysis.jl itself or develop new functionality. For that workflow:

1. Install VS Code from <https://code.visualstudio.com/Download>.
2. Install the official Julia extension.
3. Clone `https://github.com/michaelsachs/GlobalAnalysis.jl`.
4. Open the repository folder in VS Code.
5. Start a Julia REPL in VS Code.
6. Press `]` to enter the package manager and run:

```julia-repl
pkg> instantiate
```

This developer workflow is separate from the browser notebook setup above.
