```@meta
CurrentModule = MechanicalMaterialModels
```

# MechanicalMaterialModels

Implementations of mechanical (stress-strain) material models following 
the `MaterialModelsBase.jl` interface. 


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