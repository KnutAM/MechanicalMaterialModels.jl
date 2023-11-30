```@meta
CurrentModule = MechanicalMaterialModels
```

# MechanicalMaterialModels

Implementations of mechanical (stress-strain) material models following 
the `MaterialModelsBase.jl` interface. 


## Elasticity
### Linear Isotropic Elasticity
```@docs 
LinearElastic(::Val{:isotropic})
```

### Linear Anisotropic Elasticity
```@docs 
LinearElastic(::Val{:general})
LinearElastic(::Val{:cubicsymmetry})
```

## Plasticity
```@docs
Plastic
```

### Yield criteria
```@docs
VonMises
DruckerPrager
```

### Isotropic hardening
```@docs
Voce
Swift
```

### Kinematic hardening
```@docs
ArmstrongFrederick
Delobelle
OhnoWang
```

### Overstress functions
```@docs
RateIndependent
NortonOverstress
```

## Internal functions
The following functions can be relevant to use
outside the package as well, but are not part of the 
public API.
```@docs
vonmises
macaulay
```

## Developer Documentation
The following functions and their docstrings may be helpful 
for understanding the code.
```@docs
maketuple_or_nothing
get_promoted_type
baseof
vector_residual!
```