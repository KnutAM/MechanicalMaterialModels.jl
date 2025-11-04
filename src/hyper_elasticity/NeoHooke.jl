
"""
    NeoHooke(; G)

The incompressible neo-Hookean formulation with shear modulus `G` defined by the potential
```math
\\varPsi(\\boldsymbol{C}) = \\frac{G}{2} \\left[ \\frac{\\mathrm{tr}(\\boldsymbol{C})}{\\sqrt[3]{\\det(\\boldsymbol{C})}} - 3\\right]
```
where ``\\boldsymbol{C}`` is the Right Cauchy-Green deformation tensor. 

Note that sometimes, the division by ``\\sqrt[3]{\\det(\\boldsymbol{C})}`` is omitted, 
which has no influence if `C` truly is incompressible, i.e. ``\\det(\\boldsymbol{C}) = 1``. 
"""
@kwdef struct NeoHooke{T} <: AbstractHyperElastic
    G::T
end

compute_potential(m::NeoHooke, C::SymmetricTensor) = (m.G / 2) * (tr(C) / cbrt(det(C)) - 3)
MMB.get_params_eltype(::NeoHooke{T}) where {T} = T

"""
    CompressibleNeoHooke(; G, K)

A compressible neo-Hookean formulation defined by the potential
```math
\\varPsi(\\boldsymbol{C}) = 
\\frac{G}{2} \\left[ \\frac{\\mathrm{tr}(\\boldsymbol{C})}{\\sqrt[3]{\\det(\\boldsymbol{C})}} - 3\\right]
+ \\varPsi_\\mathrm{vol}\\left(\\sqrt{\\det(\\boldsymbol{C})}\\right), 
\\quad 
\\varPsi_\\mathrm{vol}(J) = \\frac{K}{2} \\left[ J - 1 \\right]^2
```

Note that there are many different variations of this model considering the volumetric part, ``\\varPsi_\\mathrm{vol}``.
"""
@kwdef struct CompressibleNeoHooke{T} <: AbstractHyperElastic
    G::T
    K::T
end

function compute_potential(m::CompressibleNeoHooke, C::SymmetricTensor)
    detC = det(C)
    ΨG = (m.G / 2) * (tr(C) / cbrt(detC) - 3)
    ΨK = (m.K / 2) * (sqrt(detC) - 1)^2
    return ΨG + ΨK
end

MMB.get_params_eltype(::CompressibleNeoHooke{T}) where {T} = T
