using Documenter
using ATTIRE

makedocs(
    sitename = "ATTIRE",
    format = Documenter.HTML(),
    modules = [ATTIRE],
    pages = [
    "Home" => "index.md",
    "Manual" => Any[
        "Guide" => "man/guide.md",
        "man/syntax.md",
    ],
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

