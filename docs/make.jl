using MechanicalMaterialModels
using Documenter

DocMeta.setdocmeta!(MechanicalMaterialModels, :DocTestSetup, :(using MechanicalMaterialModels); recursive=true)

makedocs(;
    modules=[MechanicalMaterialModels],
    authors="Knut Andreas Meyer and contributors",
    repo="https://github.com/KnutAM/MechanicalMaterialModels.jl/blob/{commit}{path}#{line}",
    sitename="MechanicalMaterialModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://KnutAM.github.io/MechanicalMaterialModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/KnutAM/MechanicalMaterialModels.jl",
    devbranch="main",
)
