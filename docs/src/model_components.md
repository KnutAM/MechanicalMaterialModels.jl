```@meta
CurrentModule = MechanicalMaterialModels
```
# Model Components

## Yield criteria
```@docs
VonMises
DruckerPrager
```

## Isotropic hardening
```@docs
Voce
Swift
```

## Kinematic hardening
```@docs
ArmstrongFrederick
```
The evolution equation according to Armstrong-Frederick hardening 
for $\boldsymbol{b}$ is given as
```math
g_{\mathrm{kin},i}(\nu, \boldsymbol{\beta}_i) = -\boldsymbol{\nu} + \frac{3\boldsymbol{\beta}_i^\mathrm{T}}{2\beta_\infty}
```
noting that for small strains, $\boldsymbol{\beta}_i = \boldsymbol{\beta}_i^\mathrm{T}$
is symmetric. 

Temp ref: doi: 10.1179/096034007X207589

```@docs
Delobelle
```
The evolution equation according to Delobelle hardening 
for $\boldsymbol{b}$ is given as
```math
g_{\mathrm{kin},i}(\nu, \boldsymbol{\beta}_i) 
= - \boldsymbol{\nu} + \delta\frac{3\boldsymbol{\beta}_i^\mathrm{T}}{2\beta_\infty}
+ \left[1 - \delta\right] \frac{\boldsymbol{\nu}:\boldsymbol{\beta}_i}{\beta_\infty}\boldsymbol{\nu}
```
noting that for small strains, $\boldsymbol{\beta}_i = \boldsymbol{\beta}_i^\mathrm{T}$
is symmetric. 

Temp ref: doi: 10.1016/S0749-6419(95)00001-1

```@docs
OhnoWang
```
The evolution equation according to Ohno-Wang hardening
for $\boldsymbol{b}$ is given as
```math
g_{\mathrm{kin},i}(\nu, \boldsymbol{\beta}_i) 
= - \boldsymbol{\nu} 
  + \frac{3\boldsymbol{\beta}_i^\mathrm{T}}{2\beta_\infty} 
    \frac{\langle \boldsymbol{\nu}:\boldsymbol{\beta}_i \rangle}{\beta_\infty}
    \left[\frac{\beta_i^\mathrm{vM}}{\beta_\infty}\right]^m
```
where $\langle x \rangle = \mathrm{max}(x, 0)$ is the Macaulay bracket and
$\beta_i^\mathrm{vM}$ is the effective von Mises backstress, see [`vonmises`](@ref).
Note that for small strains, $\boldsymbol{\beta}_i = \boldsymbol{\beta}_i^\mathrm{T}$
is symmetric. 

Temp ref: 10.1016/0749-6419(93)90042-O

## Overstress functions
```@docs
RateIndependent
NortonOverstress
```
