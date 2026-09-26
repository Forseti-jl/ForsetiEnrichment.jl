using Documenter
using ForsetiEnrichment

DocMeta.setdocmeta!(ForsetiEnrichment, :DocTestSetup, :(using ForsetiEnrichment); recursive = true)

makedocs(;
    sitename = "ForsetiEnrichment.jl",
    modules = [ForsetiEnrichment],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/Forseti-jl/ForsetiEnrichment.jl.git",
    devbranch = "main",
)
