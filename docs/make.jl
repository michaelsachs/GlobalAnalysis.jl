using Documenter
using GlobalAnalysis

makedocs(
    sitename = "GlobalAnalysis",
    format = Documenter.HTML(),
    modules = [GlobalAnalysis]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
