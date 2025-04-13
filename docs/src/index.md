```@meta
CurrentModule = MechanicalMaterialModels
```

# MechanicalMaterialModels

Implementations of mechanical (stress-strain) material models following 
the [`MaterialModelsBase.jl`](https://github.com/KnutAM/MaterialModelsBase.jl) 
interface, and to use this package that and 
[`Newton.jl`](https://github.com/KnutAM/Newton.jl) must be installed:
```julia
using Pkg
Pkg.add(url="https://github.com/KnutAM/MaterialModelsBase.jl")
Pkg.add(url="https://github.com/KnutAM/Newton.jl")
Pkg.add(url="https://github.com/KnutAM/MechanicalMaterialModels.jl")
```

See [``'s documentation](https://knutam.github.io/MaterialModelsBase.jl/dev/) for details
on how to use materials.

## Examples of available material models
* [`LinearElastic`](@ref): Isotropic and anisotropic linear elasticity for small strains
* [`GeneralizedMaxwell](@ref): The standard generalized linear viscoelasticity model with arbitrarily many maxwell-chains.
* [`Plastic`](@ref): Small strain plasticity supporting
  - Different yield surfaces, e.g. von Mises and Drucker-Prager
  - Different isotropic and kinematic hardening laws
  - Both viscoplasticity with overstress functions and standard rate-indepenent plasticity.
* [`CompressibleNeoHooke`](@ref): One of several hyperelasticity models included for finite strains
* [`FiniteStrainPlastic`](@ref): Finite strain plasticity with the same features as [`Plastic`](@ref), but for finite strains.
