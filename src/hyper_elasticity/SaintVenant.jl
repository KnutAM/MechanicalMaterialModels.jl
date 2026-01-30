
"""
    SaintVenant(elastic::LinearElastic)

The Saint-Venant formulation defined by the potential
```math
\\varPsi(\\boldsymbol{E}) = \\frac{1}{2}\\boldsymbol{E}:\\mathsf{\\boldsymbol{C}}:\\boldsymbol{E}
```
where ``\\boldsymbol{E}`` is the Green-Lagrange strain tensor and ``\\mathsf{\\boldsymbol{C}}`` is the 
elastic stiffness tensor defined in `elastic`.
"""
struct SaintVenant{E} <: AbstractHyperElastic
    elastic::E
end

function compute_potential(m::SaintVenant, C::SymmetricTensor)
    E = 0.5 * (C - one(C))
    return 1/2 * E ⊡ m.elastic.C ⊡ E
end

MMB.get_vector_eltype(m::SaintVenant) = MMB.get_vector_eltype(m.elastic)
