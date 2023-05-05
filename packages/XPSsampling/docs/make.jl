using Documenter, DocumenterTools
using XPSsampling

makedocs(
    sitename = "XPSsampling",
    format = Documenter.HTML(mathengine = MathJax3(Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
            "tags" => "ams",
            "packages" => ["base", "ams", "autoload", "physics"],
        ),
    )),
    ), #mathengine
    modules = [XPSsampling],
    pages = [
    "Home" => "index.md",
    ],
    doctest = true,
)

# deploydocs(
#     repo = "<repository url>"
# )

