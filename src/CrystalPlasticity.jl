"""
    CrystalPlasticity(;crystal, elastic, yield, q, h0, h∞, ζ, Hkin, β∞, overstress)

!!! note Experimental API
    The current implementation is one very specific model, and it is expected that it will be 
    generalized in the future to become more modular.

The parameters for this model are, 
* `crystal::Crystallography`: Which crystal structure should be used, e.g. [`BCC`](@ref) or [`FCC`](@ref).
* `elastic::LinearElastic`: Which elastic law to be used. 
* `yield::Number`: Initial yield limit in all slip systems
* `q::Number`: Cross hardening factor, q=0 is no cross hardening, q = 1 is full
* `h0::Number`: Initial hardening modulus
* `h∞::Number`: Saturated hardening modulus
* `ζ::Number`: Hardening saturation rate
* `Hkin::Number`: Kinematic hardening modulus
* `β∞::Number`: Kinematic saturation stress
* `overstress`: An overstress function (`<:Overstress`) or `RateIndependent`.

See [Meyer (2020)](https://doi.org/10.1016/j.ijsolstr.2020.04.037) for a description of this model
in a finite strain framework. 
"""
@kwdef struct CrystalPlasticity{C, E, T, OS, T_tol} <: AbstractMaterial
    crystal::C          # E.g. BCC, FCC, etc.
    elastic::E          # Elastic definition
    yield::T            # Initial yield limit
    # Isotropic hardening parameters
    q::T    # Cross hardening factor (q=0 => no cross hardening, q=1 => "full")
    h0::T   # Initial hardening modulus 
    h∞::T   # Saturated hardening modulus
    ζ::T    # Hardening saturation rate
    # Kinematic hardening parameters
    Hkin::T # Kinematic hardening modulus 
    β∞::T   # Kinematic saturation stress
    overstress::OS      # Overstress function
    maxiter::Int = 100
    tolerance::T_tol = sqrt(eps(typeof(q)))
end

struct CrystalPlasticityState{T, nslip} <: MMB.AbstractMaterialState
    γ_acc::T
    ϵp::SymmetricTensor{2, 3, T, 6}
    κ::SVector{nslip, T}
    β::SVector{nslip, T}
end

function MMB.initial_material_state(mat::CrystalPlasticity)
    T = Float64
    v = zero(SVector{get_num_slipsystems(mat.crystal), T})
    return CrystalPlasticityState(zero(T), zero(SymmetricTensor{2,3,T}), v, v)
end

function calculate_plastic_strain(mat::CrystalPlasticity, Δγ::SVector, old::CrystalPlasticityState, ν_trial)
    return old.ϵp + symmetric(sum(Δγ .* ν_trial .* get_slip_dyads(mat.crystal)))
end

function calculate_stress(mat::CrystalPlasticity, Δγ::SVector, old::CrystalPlasticityState, ϵ, ν_trial)
    ϵp = calculate_plastic_strain(mat, Δγ, old, ν_trial)
    return calculate_stress(mat.elastic, ϵ - ϵp)
end

function residual(mat::CrystalPlasticity, Δγ::SVector{nslip}, old::CrystalPlasticityState, ϵ::SymmetricTensor{2,3}, ν_trial::SVector{nslip}, Δt) where {nslip}
    σ = calculate_stress(mat, Δγ, old, ϵ, ν_trial)
    sum_Δγ = sum(Δγ)
    γ_acc = old.γ_acc + sum_Δγ
    h = mat.h∞ + (mat.h0 - mat.h∞) * exp(-mat.ζ * γ_acc)
    cross_hard = mat.q * sum_Δγ
    slip_dirs = get_slip_directions(mat.crystal)
    slip_planes = get_slip_planes(mat.crystal)
    return map(
        (dγ, κ_old, β_old, ν_tr, s, m) -> (
            τ = s ⋅ σ ⋅ m;
            κ = κ_old + h * (cross_hard + (1 - mat.q) * dγ);
            β = mat.β∞ * ν_tr + exp(-mat.Hkin * dγ / mat.β∞) * (β_old - mat.β∞ * ν_tr);
            Φ = abs(τ - β) - (mat.yield + κ);
            yield_residual(mat.overstress, Φ, dγ, Δt, mat.yield + κ)
            ),
        Δγ, old.κ, old.β, ν_trial, slip_dirs, slip_planes)
end

function calculate_state(mat::CrystalPlasticity, Δγ::SVector{nslip}, old::CrystalPlasticityState, ν_trial::SVector{nslip}) where {nslip}
    sum_Δγ = sum(Δγ)
    γ_acc = old.γ_acc + sum_Δγ
    h = mat.h∞ + (mat.h0 - mat.h∞) * exp(-mat.ζ * γ_acc)
    cross_hard = mat.q * sum_Δγ
    ϵp = calculate_plastic_strain(mat, Δγ, old, ν_trial)
    κ = map((κ_old, dγ) -> κ_old + h * (cross_hard + (1 - mat.q) * dγ), old.κ, Δγ)
    β = map((β_old, dγ, ν_tr) -> mat.β∞ * ν_tr + exp(-mat.Hkin * dγ / mat.β∞) * (β_old - mat.β∞ * ν_tr), old.β, Δγ, ν_trial)
    return CrystalPlasticityState(γ_acc, ϵp, κ, β)
end

function MMB.material_response(mat::CrystalPlasticity, ϵ::SymmetricTensor{2,3}, old::CrystalPlasticityState, Δt, _, _)
    σ_trial = calculate_stress(mat.elastic, ϵ - old.ϵp)
    τ_trial_red = map((β, s, m) -> s ⋅ σ_trial ⋅ m - β, old.β, get_slip_directions(mat.crystal), get_slip_planes(mat.crystal))
    if all(abs(τred) < (mat.yield + κ) for (τred, κ) in zip(τ_trial_red, old.κ))
        return σ_trial, mat.elastic.C, old
    else # Plastic 
        Δγ0 = zero(typeof(τ_trial_red))
        ν_trial = sign.(τ_trial_red)
        rf(x) = residual(mat, x, old, ϵ, ν_trial, Δt)
        Δγ, ∂r∂x, converged = newtonsolve(rf, Δγ0; maxiter = mat.maxiter, tol = mat.tolerance)
        if converged
            σ = calculate_stress(mat, Δγ, old, ϵ, ν_trial)
            # R(Δγ(ϵ), ϵ) = 0 => ∂R/∂x ∂γ/∂ϵ = - ∂R/∂ϵ
            # dσ/dϵ = ∂σ/∂ϵ + ∂σ/∂γ ∂γ/∂ϵ = ∂σ/∂ϵ - ∂σ/∂γ [∂R/∂x]⁻¹ ∂R/∂ϵ
            ∂σ∂γ = ForwardDiff.jacobian(γ -> tomandel(SVector, calculate_stress(mat, γ, old, ϵ, ν_trial)), Δγ)
            ∂r∂ϵ = ForwardDiff.jacobian(e -> residual(mat, Δγ, old, frommandel(SymmetricTensor{2,3}, e), ν_trial, Δt), tomandel(SVector, ϵ))
            dσdϵ = mat.elastic.C - frommandel(SymmetricTensor{4,3}, ∂σ∂γ * (∂r∂x \ ∂r∂ϵ))
            return σ, dσdϵ, calculate_state(mat, Δγ, old, ν_trial)
        else
            throw(MMB.NoLocalConvergence("$(typeof(mat)): newtonsolve! didn't converge, ϵ = ", ϵ))
        end
    end
end
