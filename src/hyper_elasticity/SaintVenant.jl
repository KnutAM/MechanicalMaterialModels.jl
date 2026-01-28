
"""
    SaintVenant(; G, K)

The-Saint Venant formulation with shear modulus `G` and bulk modulus `K` defined by the potential
```math
\\varPsi(\\boldsymbol{E}) = G\\boldsymbol{E}_{dev}:\\boldsymbol{E}_{dev} \\frac{K}{2} \\mathrm{tr}(\\boldsymbol{E})^2
```
where ``\\boldsymbol{E}`` is the Green-Lagrange strain tensor. 

"""
@kwdef struct SaintVenant{T} <: AbstractHyperElastic
    G::T
    K::T
end

function compute_potential(m::SaintVenant, C::SymmetricTensor)
    E = 0.5 * (C - one(C))
    return m.G * dev(E)⊡dev(E)+m.K/2 * (tr(E))^2
end
MMB.get_vector_eltype(::SaintVenant{T}) where {T} = T