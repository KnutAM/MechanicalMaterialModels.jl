# Material-specific implementations of `MaterialModelsBase.stress_from_state`
# (https://github.com/KnutAM/MaterialModelsBase.jl/pull/21), which upstreamed the
# generic postprocessing machinery originally prototyped in this file (see
# https://github.com/KnutAM/MechanicalMaterialModels.jl/pull/13): the generic
# `NoMaterialState` fallbacks, `FrozenStressMaterial`, the generic
# reduced-dimensional fallback, and `ReducedStressState` delegation now all live
# in `MaterialModelsBase.jl` itself. Only the material-specific methods below -
# either required (no `NoMaterialState`-based default exists) or a cheaper,
# non-autodiff alternative to the generic fallback - remain here.
#
# NOTE: while MaterialModelsBase.jl#21 is not yet merged, Project.toml and
# docs/Project.toml temporarily point MaterialModelsBase at the branch
# implementing it (knutambot/MaterialModelsBase.jl#cb/calculate_current_stress).
# Revert both `[sources]` entries (and tighten the `MaterialModelsBase` `[compat]`
# bound to the first release containing `stress_from_state`) once that PR merges
# and is released.

# LinearElastic.jl: stress-only, no gradient (material_response would compute one,
# via `m.C`, that `stress_from_state`'s generic `NoMaterialState` fallback would
# otherwise compute via `material_response` and discard).
MMB.stress_from_state(m::LinearElastic, ϵ::SymmetricTensor{2,3}, ::MMB.NoMaterialState) = calculate_stress(m, ϵ)

# HyperElastic.jl: stress-only. Computing `S = 2 ∂Ψ/∂C` (once) is unavoidable to get
# the stress at all, but the generic `NoMaterialState` fallback's `material_response`
# call additionally differentiates through that once more to get the tangent, which
# `stress_from_state` doesn't need.
MMB.stress_from_state(m::AbstractHyperElastic, F::Tensor{2,3}, ::MMB.NoMaterialState) = F ⋅ compute_stress(m, tdot(F))

# Plastic.jl
# Reduced-dimensional stress states (e.g. PlaneStress) need no dedicated method
# here: MaterialModelsBase's generic fallback autodiffs through this full-dim
# method and gives the same result (this formula is linear in ϵ, so the
# autodiff-derived tangent is exact, same as the elastic stiffness itself).
function MMB.stress_from_state(m::Plastic, ϵ::SymmetricTensor{2,3}, state::PlasticState)
    return calculate_stress(m.elastic, ϵ - state.ϵp)
end

# ViscoElastic.jl
function MMB.stress_from_state(m::GeneralizedMaxwell, ϵ::SymmetricTensor{2,3}, state::GeneralizedMaxwellState)
    σ0 = calculate_stress(m.base, ϵ)
    return mapreduce((c, ϵv) -> 2 * c.G * (dev(ϵ) - ϵv), +, m.chains, state.ϵv; init=σ0)
end

# FiniteStrainPlastic.jl
# `calculate_PKstress(m, state, F)` already computes the frozen-state (converged
# `state.Fp`, no Newton re-solve) 1st Piola-Kirchhoff stress; it is used
# internally for the elastic-predictor branch of `material_response`. No
# reduced-dimensional method is needed: MaterialModelsBase's generic fallback
# (autodiff-ing through this method via its own `FrozenStressMaterial`) covers it.
MMB.stress_from_state(m::FiniteStrainPlastic, F::Tensor{2,3}, state::FiniteStrainPlasticState) = calculate_PKstress(m, state, F)

# RotatedMaterial.jl
function _stress_from_state_rotated(rm::RotatedMaterial, ϵ::SymmetricTensor{2,3}, state)
    θ = norm(rm.rotation)
    ϵ_rot = rotate(ϵ, rm.rotation, -θ)
    σ_rot = MMB.stress_from_state(rm.material, ϵ_rot, state)
    return rotate(σ_rot, rm.rotation, θ)
end
MMB.stress_from_state(rm::RotatedMaterial, ϵ::SymmetricTensor{2,3}, state) = _stress_from_state_rotated(rm, ϵ, state)
# Disambiguates against MaterialModelsBase's own `(AbstractMaterial, ϵ,
# ::NoMaterialState)` fallback, which would otherwise be equally specific when
# `rm.material` is stateless.
MMB.stress_from_state(rm::RotatedMaterial, ϵ::SymmetricTensor{2,3}, state::MMB.NoMaterialState) = _stress_from_state_rotated(rm, ϵ, state)
