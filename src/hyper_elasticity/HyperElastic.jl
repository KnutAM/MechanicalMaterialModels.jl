abstract type AbstractHyperElastic <: AbstractMaterial end

"""
    compute_potential(m::AbstractHyperElastic, C::SymmetricTensor)

Compute the potential Ψ(C) given the Right-Cauchy-Green deformation tensor `C` 
for the hyperelastic material `m`.
"""
function compute_potential end

"""
    compute_stress(m::AbstractHyperElastic, C::SymmetricTensor)

Compute the 2nd PiolaKirchhoff stress, `S = 2 ∂Ψ/∂C`, for the potential
Ψ defined by `m` for the Right-Cauchy-Green deformation tensor `C` 
"""
compute_stress(m::AbstractHyperElastic, C::SymmetricTensor) = 2*gradient(x -> compute_potential(m, x), C)

"""
    compute_tangent(m::AbstractHyperElastic, C::SymmetricTensor)

Compute the tangent stiffness, `∂S/∂E = 4 ∂²Ψ/∂C²`, for the potential `Ψ`
defined by `m` for the Right-Cauchy-Green deformation tensor `C`.
"""
compute_tangent(m::AbstractHyperElastic, C::SymmetricTensor) = 2*gradient(x -> compute_stress(m, x), C)

"""
    compute_stress_and_tangent(m::AbstractHyperElastic, C::SymmetricTensor)

Compute both the 2nd Piola-Kirchhoff stress, `S`, and the tangent stiffness, `∂S/∂E = 4 ∂²Ψ/∂C²`, 
for the potential `Ψ` defined by `m` for the Right-Cauchy-Green deformation tensor `C`.
This is normally more efficient than computing the stress and tangents invidividually.
"""
function compute_stress_and_tangent(m::AbstractHyperElastic, C::SymmetricTensor)
    ∂S∂C, S = gradient(x -> compute_stress(m, x), C, :all)
    return S, 2*∂S∂C
end

function MMB.material_response(m::AbstractHyperElastic, F::Tensor{2,3}, args...)
    C = tdot(F)
    S, ∂S∂E = compute_stress_and_tangent(m, C)
    P = F ⋅ S 
    I2 = one(F)
    ∂P∂F = F ⋅ ∂S∂E ⊡ otimesu(F',I2) + otimesu(I2, S)
    return P, ∂P∂F, NoMaterialState()
end

MMB.get_tensorbase(::AbstractHyperElastic) = Tensor{2,3}
