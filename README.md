# MechanicalMaterialModels

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/MechanicalMaterialModels.jl/dev/)
[![Build Status](https://github.com/KnutAM/MechanicalMaterialModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/MechanicalMaterialModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/MechanicalMaterialModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/MechanicalMaterialModels.jl)

`MechanicalMaterialModels.jl` provides common material models following the
[MaterialModelsBase.jl](https://github.com/KnutAM/MaterialModelsBase.jl) interface. See the documentation for the available models. 

## Installation
The `MechanicalMaterialModels.jl` package depends on the following unregistered packages: 

* [MaterialModelsBase.jl](https://github.com/KnutAM/MaterialModelsBase.jl)
* [Newton.jl](https://github.com/KnutAM/Newton.jl)

Those packages, as well as `MechanicalMaterialModels.jl` itself, is available in [knutamregistry](https://github.com/KnutAM/knutamregistry). 
After adding this registry, the `MechanicalMaterialModels.jl` and its dependencies can be installed as other registered packages, i.e.
```julia
using Pkg
Pkg.add("MechanicalMaterialModels")
```

If the registry is not added, it is necessary to install all 
non-registered required packages, i.e.
```julia
using Pkg
Pkg.add(url="https://github.com/KnutAM/Newton.jl")
Pkg.add(url="https://github.com/KnutAM/MaterialModelsBase.jl")
Pkg.add(url="https://github.com/KnutAM/MechanicalMaterialModels.jl")
```
