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

A material-model developer only needs to implement the full-dimensional method,
`calculate_current_stress(m::MyMaterial, ϵ, state::MyMaterialState)`. Support for
a reduced-dimensional stress state (via `ReducedStressState`) then follows
automatically from a generic fallback, which rides `MaterialModelsBase`'s
existing stress-state Newton iteration (e.g. `PlaneStress`) using an internal
`FrozenStressMaterial` wrapper, with the tangent obtained by automatic
differentiation. A specific reduced-dimensional method only needs to be added
when a cheaper, non-autodiff alternative exists (as done here for
`Plastic` and for stateless/`NoMaterialState` materials).

Currently supported materials: [`LinearElastic`](@ref), [`NeoHooke`](@ref),
[`CompressibleNeoHooke`](@ref), and [`SaintVenant`](@ref) (each with a
dedicated, gradient-free implementation), [`Plastic`](@ref),
[`FiniteStrainPlastic`](@ref), [`GeneralizedMaxwell`](@ref), and
[`RotatedMaterial`](@ref) wrapping any of these. Reduced-dimensional support
(via `ReducedStressState`) works for all of the above, generically for
`GeneralizedMaxwell` and for `RotatedMaterial` wrapping a small-strain
material (via the generic fallback), and with a dedicated
non-autodiff implementation for stateless materials and for `Plastic`.

!!! note "Not (yet) supported"
    `CrystalPlasticity` (small-strain, despite referencing a finite-strain
    framework in its docstring) has no `calculate_current_stress` method at
    all yet. `RotatedMaterial` wrapping a finite-strain material already
    errors in `RotatedMaterial`'s own `material_response` (a hard
    `::SymmetricTensor{2,3}` type assertion), independently of this function.
"""
function calculate_current_stress end

# LinearElastic.jl: stress-only, no gradient (material_response would compute one,
# via `m.C`, that `calculate_current_stress` doesn't need).
calculate_current_stress(m::LinearElastic, ϵ::SymmetricTensor{2,3}, ::MMB.NoMaterialState) = calculate_stress(m, ϵ)

# HyperElastic.jl: stress-only. Computing `S = 2 ∂Ψ/∂C` (once) is unavoidable to get
# the stress at all, but `material_response` additionally differentiates through
# that once more to get the tangent, which `calculate_current_stress` doesn't need.
calculate_current_stress(m::AbstractHyperElastic, F::Tensor{2,3}, ::MMB.NoMaterialState) = F ⋅ compute_stress(m, tdot(F))

# Plastic.jl
function calculate_current_stress(m::Plastic, ϵ::SymmetricTensor{2,3}, state::PlasticState)
    return calculate_stress(m.elastic, ϵ - state.ϵp)
end

function calculate_current_stress(stress_state::MMB.AbstractStressState, m::Plastic, ϵ, state::PlasticState)
    # Expand the (possibly reduced) total strain to 3d before removing the plastic
    # strain: for non-iterative states (e.g. PlaneStrain) the zero-padded
    # out-of-plane *total* strain is exact by definition of the state, whereas
    # reducing `state.ϵp` first would incorrectly discard its out-of-plane part.
    # This avoids autodiff entirely, by delegating to `m.elastic`'s own analytic
    # stress-state response.
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

# Wraps a frozen-state stress formula, `f`, mapping a strain (`SecondOrderTensor{3}`,
# i.e. `Tensor{2,3}` or `SymmetricTensor{2,3}`) to a stress (at fixed history/internal
# variables) as an `AbstractMaterial`, so that it can ride MaterialModelsBase's
# existing stress-state Newton iteration (e.g. for `PlaneStress`). The tangent needed
# for that iteration is obtained via automatic differentiation. This is what powers
# the generic reduced-dimensional fallback of `calculate_current_stress` below.
struct FrozenStressMaterial{F} <: AbstractMaterial
    f::F
end
function MMB.material_response(fm::FrozenStressMaterial, strain::SecondOrderTensor{3}, old::MMB.AbstractMaterialState, args::Vararg{Any,N}) where {N}
    dσdϵ, σ = Tensors.gradient(fm.f, strain, :all)
    return σ, dσdϵ, old
end

# Generic reduced-dimensional fallback: as long as `calculate_current_stress(m, ϵ,
# state)` (full-dimensional) is implemented for `m`, this makes `ReducedStressState`
# support "just work", by autodiff-ing through it. More specific methods above/below
# (e.g. for `Plastic` or `NoMaterialState`) take precedence when a cheaper,
# non-autodiff alternative exists.
function calculate_current_stress(stress_state::MMB.AbstractStressState, m::AbstractMaterial, strain, state::MMB.AbstractMaterialState)
    frozen = FrozenStressMaterial(e -> calculate_current_stress(m, e, state))
    σ, _, _, _ = MMB.material_response(stress_state, frozen, strain, MMB.NoMaterialState{eltype(strain)}())
    return σ
end

# Reduced-dimensional fast path for stateless materials: avoids the autodiff in the
# generic fallback above by delegating directly to `material_response`'s own
# (analytic, for `LinearElastic`) stress-state handling.
function calculate_current_stress(stress_state::MMB.AbstractStressState, m::AbstractMaterial, strain, state::MMB.NoMaterialState)
    σ, _, _, _ = MMB.material_response(stress_state, m, strain, state)
    return σ
end

# RotatedMaterial.jl
function calculate_current_stress(rm::RotatedMaterial, ϵ::SymmetricTensor{2,3}, state)
    θ = norm(rm.rotation)
    ϵ_rot = rotate(ϵ, rm.rotation, -θ)
    σ_rot = calculate_current_stress(rm.material, ϵ_rot, state)
    return rotate(σ_rot, rm.rotation, θ)
end

# ReducedStressState (MaterialModelsBase.jl)
function calculate_current_stress(rss::MMB.ReducedStressState, ϵ, state)
    return calculate_current_stress(rss.stress_state, rss.material, ϵ, state)
end
