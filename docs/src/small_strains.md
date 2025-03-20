```@meta
CurrentModule = MechanicalMaterialModels
```
# Models for Small Strains
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

## [Plasticity](@id small_strain_plasticity)
```@docs
Plastic
```

### [Theory](@id ssp_theory)
While the exact model response is given by the laws in `elastic`, `yield`, `isotropic`, `kinematic`, and `overstress`,
the generic model equations are described below. 

The stress is calculated from the elastic strains, ``\boldsymbol{\epsilon}_\mathrm{e}``, obtained via the 
additive decomposition, ``\boldsymbol{\epsilon} = \boldsymbol{\epsilon}_\mathrm{e} + \boldsymbol{\epsilon}_\mathrm{p}``. 
The elastic law is specified by `m.elastic` and is evaluated by giving it the elastic strain. 

#### [Yield Criterion](@id ssp_yield)
A yield criterion of the type 
```math
\varPhi = f\left( \boldsymbol{\sigma} - \boldsymbol{\beta} \right) - \left[Y_0 + \kappa\right]
```
is assumed. Here, ``\boldsymbol{\beta} = \sum_{i=1}^{N_\mathrm{kin}} \boldsymbol{\beta}_i`` is the total back-stress, 
and ``\kappa = \sum_{i=1}^{N_\mathrm{iso}} \kappa_i`` is the total isotropic hardening stress. The initial yield limit 
is passed to the yield criterion along with potentially other parameters. 
The evolution laws for ``\boldsymbol{\beta}_i`` and ``\kappa_i`` are given by the kinematic and isotropic hardening laws.

Associative plastic flow is used to obtain the plastic strains,
```math
\dot{\epsilon}_{\mathrm{p}} = \dot{\lambda} \left.\frac{\partial f}{\partial \boldsymbol{\sigma}}\right\vert_{\left( \boldsymbol{\sigma} - \boldsymbol{\beta} \right)}
= \dot{\lambda} \boldsymbol{\nu}
```

#### [Isotropic Hardening](@id ssp_isohard)
The isotropic hardening stress is $\kappa_i = -H_{\mathrm{iso},i} k_i$, where the evolution law for $k_i$ is
```math
\dot{k}_i = \dot{\lambda} g_{\mathrm{iso},i}(\kappa_i)
```
where $g_{\mathrm{iso},i}$ is specified by an isotropic hardening law, see [Isotropic hardening](@ref). The hardening modulus, $H_{\mathrm{iso},i}$, is a parameter given to the isotropic hardening law. 

#### [Kinematic Hardening](@id ssp_kinhard)
The back-stress is $\boldsymbol{\beta}_i = - [2/3]H_{\mathrm{kin},i} \boldsymbol{b}_i$, where the evolution law for $\boldsymbol{b}_i$ is
```math
\dot{\boldsymbol{b}}_i = \dot{\lambda} \boldsymbol{g}_{\mathrm{kin},i}(\boldsymbol{\nu}, \boldsymbol{\beta}_i)
```
where $g_{\mathrm{kin},i}$ is specified by a kinematic hardening law, see  [Kinematic hardening](@ref). The hardening modulus, $H_{\mathrm{kin},i}$, is a parameter given to the kinematic hardening law. 

#### [Rate Dependence](@id ssp_ratedep)
If `overstress=RateIndependent()`, the plastic multiplier, ``\lambda``, is obtained via the KKT-conditions,
```math
\dot{\lambda} \geq 0, \quad \varPhi \leq 0, \quad \dot{\lambda}\varPhi = 0
```
Otherwise, the overstress function, ``\eta(\varPhi)``, determines the evolution of ``\lambda`` as 
```math
\dot{\lambda} = \eta(\varPhi, (Y_0 + \kappa))
```

## [Viscoelasticity](@id small_strain_viscoelasticity)
```@docs
GeneralizedMaxwell
Maxwell
```
