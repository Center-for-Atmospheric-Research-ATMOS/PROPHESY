using Documenter, DocumenterTools
using XPSinv

makedocs(
    sitename = "XPSinv",
    format = Documenter.HTML(), #mathengine
    modules = [XPSinv],
    pages = [
    "Home" => "index.md",
    ],
    doctest = true,
)

# deploydocs(
#     repo = "<repository url>"
# )

# mathengine = MathJax3(Dict(
#         :loader => Dict("load" => ["[tex]/physics"]),
#         :tex => Dict(
#             "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
#             "tags" => "ams",
#             "packages" => ["base", "ams", "autoload", "physics"],
#         ),
#     )),