"""
    RotatedMaterial(m::AbstractMaterial, r::Vec) <: AbstractMaterial

A rotated material is wrapper around `m`, such that the material response is evaluated in a
local coordinate system defined by the rotation vector, `r`. The direction of `r` is the rotation
axis and `θ = |r|` (radians) is the rotation angle around that axis, following the standard
axis-angle (exponential map) convention. (Note that this is *not* the same as the Gibbs/Rodrigues
vector `tan(θ/2) * axis` found in some other references, despite both constructions building on
Rodrigues' rotation formula.) Specifically, the input strain is first rotated by `θ = -|r|` around
`r` before calling `material_response` with `m`. The resulting stress and stiffness are rotated
back, `θ = +|r|` around `r`, before they are returned.

This is equivalent to rotating the material parameter tensors (e.g. stiffness tensor `θ = +|r|` around `r`).

A zero rotation vector, `r = zero(Vec{3})`, is well-defined and gives the identity rotation
(i.e. `m` is evaluated directly, unrotated).

!!! note State variables are in the local coordinate system
    The state variables are not modified, and are hence defined in the local coordinate system.
    Care must therefore be taken to rotate the strain to the local coordinates when evaluating the
    responses, and the output of those evaluations should be rotated back to the global coordinates.

"""
struct RotatedMaterial{M <: AbstractMaterial, RT <: Vec} <: AbstractMaterial
    material::M
    rotation::RT
end

for op in ( :initial_material_state, :allocate_material_cache,
            :get_tensorbase, :get_num_statevars,
            :allocate_differentiation_output)
    @eval @inline MMB.$op(rm::RotatedMaterial) = MMB.$op(rm.material)
end

# Rodrigues' rotation formula evaluated directly from the rotation vector `r`
# (`θ = |r|`, axis `= r/θ`), without ever dividing by `θ` or `norm(r)`. The
# coefficients `sin(θ)/θ` and `(1-cos(θ))/θ²` are even, analytic functions of
# `θ`, so they are smooth functions of `θ² = r⋅r` with no singularity at
# `r = 0` (unlike `Tensors.rotation_tensor`, which normalizes the axis via
# `u/norm(u)` and is therefore 0/0 at `r = 0`).
function _rotation_tensor(r::Vec{3, T}) where T
    θsq = r ⋅ r
    a, b = _rodrigues_coefficients(θsq)
    z = zero(T)
    Ω = Tensor{2, 3}((z, r[3], -r[2], -r[3], z, r[1], r[2], -r[1], z))
    return one(Ω) + a * Ω + b * (Ω ⋅ Ω)
end

function _rodrigues_coefficients(θsq)
    if θsq < 1.0e-8
        a = 1 - θsq / 6 * (1 - θsq / 20 * (1 - θsq / 42))
        b = (1 - θsq / 12 * (1 - θsq / 30 * (1 - θsq / 56))) / 2
    else
        θ = sqrt(θsq)
        sinθ, cosθ = sincos(θ)
        a = sinθ / θ
        b = (1 - cosθ) / θsq
    end
    return a, b
end

_rotate_by_vector(x::SymmetricTensor{2, 3}, r::Vec{3}) = symmetric(_rotate_by_vector(Tensor{2,3}(x), r))
function _rotate_by_vector(x::Tensor{2, 3}, r::Vec{3})
    R = _rotation_tensor(r)
    return R ⋅ x ⋅ R'
end
_rotate_by_vector(x::SymmetricTensor{4, 3}, r::Vec{3}) = symmetric(_rotate_by_vector(Tensor{4,3}(x), r))
function _rotate_by_vector(x::Tensor{4, 3}, r::Vec{3})
    R = _rotation_tensor(r)
    return otimesu(R, R) ⊡ x ⊡ otimesu(R', R')
end

function MMB.material_response(rm::RotatedMaterial, strain::AbstractTensor, args::Vararg{Any, N}) where {N}
    strain_rot = _rotate_by_vector(strain, -rm.rotation)
    stress_rot, stiff_rot, state = MMB.material_response(rm.material, strain_rot, args...)
    stress = _rotate_by_vector(stress_rot, rm.rotation)
    stiff = _rotate_by_vector(stiff_rot, rm.rotation)
    return stress, stiff, state
end

function MMB.stress_from_state(rm::RotatedMaterial, ϵ::SecondOrderTensor{3}, state)
    ϵ_rot = _rotate_by_vector(ϵ, -rm.rotation)
    σ_rot = MMB.stress_from_state(rm.material, ϵ_rot, state)
    return _rotate_by_vector(σ_rot, rm.rotation)
end
