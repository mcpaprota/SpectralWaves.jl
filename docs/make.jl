using SpectralWaves
using Documenter

DocMeta.setdocmeta!(SpectralWaves, :DocTestSetup, :(using SpectralWaves); recursive=true)

makedocs(;
    modules=[SpectralWaves],
    authors="Maciej Paprota <mapap@ibwpan.gda.pl> and contributors",
    sitename="SpectralWaves.jl",
    format=Documenter.HTML(;
        canonical="https://mcpaprota.github.io/SpectralWaves.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => "guide.md",
        "Examples" => "examples.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
)

deploydocs(;
    repo="github.com/mcpaprota/SpectralWaves.jl",
    devbranch="main",
)
