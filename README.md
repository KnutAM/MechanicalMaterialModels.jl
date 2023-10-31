# MechanicalMaterialModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KnutAM.github.io/MechanicalMaterialModels.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/MechanicalMaterialModels.jl/dev/)
[![Build Status](https://github.com/KnutAM/MechanicalMaterialModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/MechanicalMaterialModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/MechanicalMaterialModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/MechanicalMaterialModels.jl)

The `MechanicalMaterialModels.jl` package depends on the following unregistered packages: 

* [MaterialModelsBase.jl](https://github.com/KnutAM/MaterialModelsBase.jl)
* [Newton.jl](https://github.com/KnutAM/Newton.jl)

Those packages, as well as `MechanicalMaterialModels.jl` itself, is available in [knutamregistry](https://github.com/KnutAM/knutamregistry). 
After adding this registry, the `MechanicalMaterialModels.jl` and its dependencies can be installed as 
```julia
using Pkg
Pkg.add("MechanicalMaterialModels")
```

If the registry is not added, it is possible to install using 
```julia
using Pkg
Pkg.add(url="https://github.com/KnutAM/MechanicalMaterialModels.jl")
```
But then it is first required to add all packages not in the general registry manually in the same fashion. 
