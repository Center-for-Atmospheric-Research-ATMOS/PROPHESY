using Documenter
using ATTIRE

makedocs(
    sitename = "ATTIRE",
    format = Documenter.HTML(),
    modules = [ATTIRE]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
