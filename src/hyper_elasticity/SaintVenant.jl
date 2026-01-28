
"""
    SaintVenant(; G, K)

The-Saint Venant formulation with shear modulus `G` and bulk modulus `K` defined by the potential
```math
\\varPsi(\\boldsymbol{E}) = G\\boldsymbol{E}_{dev}:\\boldsymbol{E}_{dev} \\frac{K}{2} \\mathrm{tr}(\\boldsymbol{E})^2
```
where ``\\boldsymbol{E}`` is the Green-Lagrange strain tensor. 

"""
@kwdef struct SaintVenant{E} <: AbstractHyperElastic
    elastic::E
end

function compute_potential(m::SaintVenant, C::SymmetricTensor)
    E = 0.5 * (C - one(C))
    return 1/2 * E ⊡ m.elastic.C ⊡ E
end

MMB.get_vector_eltype(::SaintVenant{T}) where {T} = T