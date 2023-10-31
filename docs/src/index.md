```@meta
CurrentModule = MechanicalMaterialModels
```

# MechanicalMaterialModels

Documentation for [MechanicalMaterialModels](https://github.com/KnutAM/MechanicalMaterialModels.jl).


## Elasticity
```@docs 
LinearElastic
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
```@docs
vonmises
macaulay
```