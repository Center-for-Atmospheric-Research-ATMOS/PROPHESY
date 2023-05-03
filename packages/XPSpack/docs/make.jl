using Documenter, DocumenterTools
using XPSpack

makedocs(
    sitename = "XPSpack",
    format = Documenter.HTML(),
    modules = [XPSpack],
    pages = [
    "Home" => "index.md",
    ],
    doctest = true,
)

# deploydocs(
#     repo = "<repository url>"
# )