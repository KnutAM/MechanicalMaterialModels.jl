```@meta
CurrentModule = MechanicalMaterialModels
```
# Models for Finite Strains
Implementations of mechanical (stress-strain) material models following 
the `MaterialModelsBase.jl` interface. 


## Hyperelasticity
### Neo-Hookean Materials
```@docs 
NeoHooke
CompressibleNeoHooke
```

## [Plasticity](@id finite_strain_plasticity)
```@docs
FiniteStrainPlastic
```

### [Theory](@id fsp_theory)
The overall modeling framework and general equations for `FiniteStrainPlastic` is described below.
The exact model is specified by the provided laws as keyword arguments during construction.

The model is based on the multiplicative split of the deformation gradient, $\boldsymbol{F}$, 
into elastic, $\boldsymbol{F}_\mathrm{e}$, and plastic, $\boldsymbol{F}_\mathrm{p}$, parts,
```math
\boldsymbol{F}=\boldsymbol{F}_\mathrm{e} \boldsymbol{F}_\mathrm{p}
```
#### [Elastic law](@id fsp_elasticity)
The 2nd Piola-Kirchhoff stress, $\boldsymbol{S}_\mathrm{e}$, 
on the intermediate (elastic) configuration is then 
calculated from the hyperelastic potential, $\varPsi(\boldsymbol{C}_\mathrm{e})$,
```math
\boldsymbol{S}_\mathrm{e} = 2\frac{\partial \varPsi}{\partial \boldsymbol{C}_\mathrm{e}}
```
where the elastic Right Cauchy-Green deformation tensor, $\boldsymbol{C}_\mathrm{e} := \boldsymbol{F}_\mathrm{e}^\mathrm{T}\boldsymbol{F}_\mathrm{e}$. We further introduce the 
Mandel stress as
```math
\boldsymbol{M} = 2 \boldsymbol{C}_\mathrm{e} \frac{\partial \varPsi}{\partial \boldsymbol{C}_\mathrm{e}}
```
where $\boldsymbol{M} = \boldsymbol{F}_\mathrm{e}^\mathrm{T} \boldsymbol{\tau} \boldsymbol{F}_\mathrm{e}^\mathrm{-T}$ and $\boldsymbol{\tau}$ is the Kirchhoff stress. 

In summary, the elastic response is determined by the potential $\varPsi$, which is specified by the supplied elastic law.

#### [Yield Criterion](@id fsp_yield)
A yield criterion of the type 
```math
\varPhi = f\left( \boldsymbol{M} - \boldsymbol{M}_\mathrm{k} \right) - \left[Y_0 + \kappa\right]
```
is formulated in terms of Mandel stresses, where 
$\boldsymbol{M}_\mathrm{k}$ is the total back-stress,
and $\kappa$ the total isotropic hardening stress. 
```math
\boldsymbol{M}_\mathrm{k} = \sum_{i=1}^{N_\mathrm{kin}} \boldsymbol{M}_{\mathrm{k},i},
\quad
\kappa = \sum_{i=1}^{N_\mathrm{iso}} \kappa_i
```
How the back-stresses, $\boldsymbol{M}_{\mathrm{k},i}$, and isotropic hardening stresses, $\kappa_i$, are calculated is further described below.

The initial yield limit, $Y_0$, is passed to the yield criterion along with potentially other parameters affecting the function $f$. 

In this framework, an associative plastic flow is used to obtain the plastic deformations,
specifically
```math
\boldsymbol{L}_\mathrm{p} := \dot{\boldsymbol{F}}_\mathrm{p} \boldsymbol{F}_\mathrm{p}^{-1} = \dot{\lambda} \mathrm{\nu}, \quad \mathrm{\nu} := \frac{\partial \varPhi}{\partial \boldsymbol{M}}
```

#### [Isotropic Hardening](@id fsp_isohard)
The isotropic hardening is formulated as $\kappa_i = g_{\mathrm{iso},i}(\lambda)$,
determined by the specific isotropic hardening laws, see [Isotropic hardening](@ref).

#### [Kinematic Hardening](@id fsp_kinhard)
The back-stresses, $\boldsymbol{M}_{\mathrm{k},i}$, are calculated as
```math
\boldsymbol{M}_{\mathrm{k},i} = 2\boldsymbol{c}_{\mathrm{k},i} \frac{\partial\varPsi_{\mathrm{k},i}}{\partial \boldsymbol{c}_{\mathrm{k},i}}
```
where the deformation tensor associated with kinematic hardening,  $\boldsymbol{c}_{\mathrm{k},i} := \boldsymbol{F}_{\mathrm{k},i}^{-\mathrm{T}}\boldsymbol{F}_{\mathrm{k},i}^{-1}$. The evolution of the kinematic hardening deformation gradient, $\boldsymbol{F}_{\mathrm{k},i}$, is then given by the specific kinematic hardening law,
```math
\boldsymbol{L}_\mathrm{p} := \dot{\boldsymbol{F}}_\mathrm{p} \boldsymbol{F}_\mathrm{p}^{-1} = \dot{\lambda} \boldsymbol{g}_{\mathrm{k},i}(\boldsymbol{\nu}, \boldsymbol{M}_{\mathrm{k},i})
```
where $\boldsymbol{g}_{\mathrm{k},i}$ are determined by the specific 
hardening laws, see [Kinematic hardening](@ref).

#### [Rate Dependence](@id fsp_ratedep)
If `overstress=RateIndependent()`, the plastic multiplier, ``\lambda``, is obtained via the KKT-conditions,
```math
\dot{\lambda} \geq 0, \quad \varPhi \leq 0, \quad \dot{\lambda}\varPhi = 0
```
Otherwise, the overstress function, ``\eta(\varPhi)``, determines the evolution of ``\lambda`` as 
```math
\dot{\lambda} = \eta(\varPhi, Y_0 + \kappa)
```
