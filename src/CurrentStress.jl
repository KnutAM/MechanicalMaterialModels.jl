"""
    calculate_current_stress(m::AbstractMaterial, ϵ, state::AbstractMaterialState)
    calculate_current_stress(rss::ReducedStressState, ϵ, state::AbstractMaterialState)

Calculate the stress that is energy-conjugated to `ϵ`, consistent with the *given*
`state`, without invoking any local iteration that would advance history/internal
variables. `state` is normally the already-converged state obtained from a previous
call to `material_response` (e.g. during postprocessing, where `ϵ` may differ
slightly from the strain that produced `state`, such as an interpolated quadrature
point value).

This is a prototype for [MaterialModelsBase.jl#12](https://github.com/KnutAM/MaterialModelsBase.jl/issues/12),
exploring how such an interface would work for different material models,
including support for a reduced-dimensional stress state
(see `MaterialModelsBase.ReducedStressState`) when only the
reduced-dimensional strain is supplied. This replaces the ad hoc, per-material
`calculate_stress` dispatch previously used for postprocessing in
[FerriteAssembly.jl#94](https://github.com/KnutAM/FerriteAssembly.jl/pull/94).

Currently supported materials: [`LinearElastic`](@ref), [`NeoHooke`](@ref),
[`CompressibleNeoHooke`](@ref), and [`SaintVenant`](@ref) (all via the generic
`NoMaterialState` fallback), [`Plastic`](@ref), [`FiniteStrainPlastic`](@ref),
[`GeneralizedMaxwell`](@ref), and [`RotatedMaterial`](@ref) wrapping any of
the small-strain materials above. Reduced-dimensional support (via
`ReducedStressState`) is currently implemented for `LinearElastic`,
`NeoHooke`, `CompressibleNeoHooke`, `SaintVenant`, `Plastic`, and
`FiniteStrainPlastic`.

!!! note "Not (yet) supported"
    `CrystalPlasticity` (small-strain, despite referencing a finite-strain
    framework in its docstring), `GeneralizedMaxwell`/`RotatedMaterial` under
    `ReducedStressState`, and `RotatedMaterial` wrapping a finite-strain
    material (the latter already errors in `RotatedMaterial`'s own
    `material_response`, independently of this function).
"""
function calculate_current_stress end

# Generic fallback for genuinely stateless materials: since there is no history
# to accidentally advance, delegating to `material_response` is safe and exact.
function calculate_current_stress(m::AbstractMaterial, ϵ, state::MMB.NoMaterialState)
    σ, _, _ = MMB.material_response(m, ϵ, state)
    return σ
end

function calculate_current_stress(stress_state::MMB.AbstractStressState, m::AbstractMaterial, ϵ, state::MMB.NoMaterialState)
    σ, _, _, _ = MMB.material_response(stress_state, m, ϵ, state)
    return σ
end

# Plastic.jl
function calculate_current_stress(m::Plastic, ϵ::SymmetricTensor{2,3}, state::PlasticState)
    return calculate_stress(m.elastic, ϵ - state.ϵp)
end

function calculate_current_stress(stress_state::MMB.AbstractStressState, m::Plastic, ϵ, state::PlasticState)
    # Expand the (possibly reduced) total strain to 3d before removing the plastic
    # strain: for non-iterative states (e.g. PlaneStrain) the zero-padded
    # out-of-plane *total* strain is exact by definition of the state, whereas
    # reducing `state.ϵp` first would incorrectly discard its out-of-plane part.
    ϵ_3d = MMB.expand_tensordim(stress_state, ϵ)
    ϵₑ = ϵ_3d - state.ϵp
    σ, _, _, _ = MMB.material_response(stress_state, m.elastic, ϵₑ, MMB.initial_material_state(m.elastic))
    return σ
end

# ViscoElastic.jl
function calculate_current_stress(m::GeneralizedMaxwell, ϵ::SymmetricTensor{2,3}, state::GeneralizedMaxwellState)
    σ0 = calculate_stress(m.base, ϵ)
    return mapreduce((c, ϵv) -> 2 * c.G * (dev(ϵ) - ϵv), +, m.chains, state.ϵv; init=σ0)
end

# FiniteStrainPlastic.jl
# `calculate_PKstress(m, state, F)` already computes the frozen-state (converged
# `state.Fp`, no Newton re-solve) 1st Piola-Kirchhoff stress; it is used
# internally for the elastic-predictor branch of `material_response`.
calculate_current_stress(m::FiniteStrainPlastic, F::Tensor{2,3}, state::FiniteStrainPlasticState) = calculate_PKstress(m, state, F)

# Wraps a frozen-state stress formula (strain -> stress, at fixed history
# variables) as an `AbstractMaterial`, so that it can ride MaterialModelsBase's
# existing stress-state Newton iteration (e.g. for `PlaneStress`). The tangent
# needed for that iteration is obtained via automatic differentiation, exactly
# analogous to how `compute_stress_and_tangent` differentiates through
# `compute_stress` in `HyperElastic.jl`.
struct FrozenStressMaterial{F} <: AbstractMaterial
    f::F  # F::Tensor{2,3} -> P::Tensor{2,3}
end
function MMB.material_response(fm::FrozenStressMaterial, F::Tensor{2,3}, old::MMB.AbstractMaterialState, args...)
    dPdF, P = Tensors.gradient(fm.f, F, :all)
    return P, dPdF, old
end

function calculate_current_stress(stress_state::MMB.AbstractStressState, m::FiniteStrainPlastic, F, state::FiniteStrainPlasticState)
    frozen = FrozenStressMaterial(F_ -> calculate_PKstress(m, state, F_))
    σ, _, _, _ = MMB.material_response(stress_state, frozen, F, MMB.NoMaterialState{eltype(F)}())
    return σ
end

# RotatedMaterial.jl
function _calculate_current_stress_rotated(rm::RotatedMaterial, ϵ::SymmetricTensor{2,3}, state)
    θ = norm(rm.rotation)
    ϵ_rot = rotate(ϵ, rm.rotation, -θ)
    σ_rot = calculate_current_stress(rm.material, ϵ_rot, state)
    return rotate(σ_rot, rm.rotation, θ)
end
calculate_current_stress(rm::RotatedMaterial, ϵ::SymmetricTensor{2,3}, state) = _calculate_current_stress_rotated(rm, ϵ, state)
# Disambiguates against the `(AbstractMaterial, ϵ, ::NoMaterialState)` fallback above,
# which would otherwise be equally specific when `rm.material` is stateless.
calculate_current_stress(rm::RotatedMaterial, ϵ::SymmetricTensor{2,3}, state::MMB.NoMaterialState) = _calculate_current_stress_rotated(rm, ϵ, state)

# ReducedStressState (MaterialModelsBase.jl)
function _calculate_current_stress_reduced(rss::MMB.ReducedStressState, ϵ, state)
    return calculate_current_stress(rss.stress_state, rss.material, ϵ, state)
end
calculate_current_stress(rss::MMB.ReducedStressState, ϵ, state) = _calculate_current_stress_reduced(rss, ϵ, state)
# Disambiguates against the `(AbstractMaterial, ϵ, ::NoMaterialState)` fallback above,
# which would otherwise be equally specific when `rss.material` is stateless.
calculate_current_stress(rss::MMB.ReducedStressState, ϵ, state::MMB.NoMaterialState) = _calculate_current_stress_reduced(rss, ϵ, state)
