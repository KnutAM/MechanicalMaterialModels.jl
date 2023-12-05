"""
    FiniteStrainPlastic(;elastic, yield, isotropic, kinematic, overstress)

A finite-strain plasticity model with modular elastic laws, yield criteria, 
multiple isotropic and kinematic hardening contributions,
and either rate-independent or viscoplastic response.
# Keyword arguments
- `elastic::AbstractMaterial`\\
  Hyperelastic law, see e.g. [`CompressibleNeoHooke`](@ref)
- `yield::YieldCriterion`\\
  Yield criterion, including the initial yield limit. If `yield::Real` is given, `VonMises(yield)` is used. 
- `isotropic::Union{AbstractIsotropicHardening,Tuple}`\\
  Isotropic hardening laws, see e.g. [`Voce`](@ref)
- `kinematic::Union{AbstractKinematicHardening,Tuple}`\\
  Kinematic hardening laws, see e.g. [`ArmstrongFrederick`](@ref)
- `overstress::Union{RateIndependent,OverstressFunction}`\\
  Rate dependence, see e.g. [`NortonOverstress`](@ref)\\
  Defaults to `RateIndependent()`

# Example
```julia
m = Plastic(elastic = CompressibleNeoHooke(G=80.e3, K=160.e3),
            yield = 300.0,
            isotropic = (Voce(Hiso=-100.e3, κ∞=-100.0),Voce(Hiso=10.e3, κ∞=200.0)),
            kinematic = (ArmstrongFrederick(Hkin=200.e3, β∞=200.0),
                         OhnoWang(Hkin=1000.e3, β∞=200.0, m=3.0)),
            overstress = NortonOverstress(;tstar=1.0, nexp=2.0))
```

"""
struct FiniteStrainPlastic{E,YLD,IsoH,KinH,OS} <: AbstractMaterial
    elastic::E          # Elastic definition
    yield::YLD          # Initial yield limit
    isotropic::IsoH     # Tuple of isotropic hardening definitions
    kinematic::KinH     # Tuple of kinematic hardening definitions
    overstress::OS      # Overstress function
end
function FiniteStrainPlastic(;elastic::AbstractHyperElastic, yield, isotropic=nothing, kinematic=nothing, overstress=RateIndependent())
    iso = maketuple_or_nothing(isotropic)
    kin = maketuple_or_nothing(kinematic)
    yld = default_yield_criteria(yield)
    return FiniteStrainPlastic(elastic, yld, iso, kin, overstress)
end


get_base_numbertype(m::FiniteStrainPlastic) = typeof(initial_yield_limit(m.yield))

# Definition of material state
struct FiniteStrainPlasticState{NKin,NIso,T} <: AbstractMaterialState
    Fp::Tensor{2,3,T,9}
    κ::NTuple{NIso, T}
    Fk::NTuple{NKin, Tensor{2,3,T,9}}
end

function MMB.initial_material_state(m::FiniteStrainPlastic)
    T = get_base_numbertype(m)
    I2 = one(Tensor{2,3,T})
    FiniteStrainPlasticState(
        I2,                                  # Fp
        map(get_initial_value, m.isotropic), # κ
        map(Returns(I2), m.kinematic))       # Fk
end

#=
# Mandel conversion of state only required for differentiation, to be done later.
Tensors.n_components(::Type{<:PlasticState{NK, NI}}) where {NK, NI} = NI + 9*(1+NK)

MMB.get_num_statevars(m::Plastic) = length(m.isotropic) + 6*(1 + length(m.kinematic))

function Tensors.tomandel!(v::AbstractVector, s::PlasticState{NKin,NIso,Tϵₚ,Tκ,Tβ}; offset=0) where {NKin,NIso,Tϵₚ,Tκ,Tβ}
    tomandel!(v, s.ϵₚ; offset=offset)
    v[(offset+7):(offset+6+NIso)] .= s.κ
    for i=1:NKin
        tomandel!(v, s.β[i], offset=offset+NIso+6*i) # 6 components in SymmetricTensor{2,3}
    end
    return v
end
function Tensors.tomandel(s::PlasticState{NKin,NIso,Tϵₚ,Tκ,Tβ}) where {NKin,NIso,Tϵₚ,Tκ,Tβ}
    v = zeros(promote_type(Tϵₚ,Tκ,Tβ), Tensors.n_components(PlasticState{NKin,NIso}))
    return tomandel!(v, s)
end

function Tensors.frommandel(::Type{<:PlasticState{NKin,NIso}}, v::AbstractVector{Tv}; offset=0) where {Tv,NKin,NIso}
    ϵₚ = frommandel(SymmetricTensor{2,3}, v; offset=offset)
    κ = ntuple(i->v[offset+6+i], NIso)
    β = ntuple(i->frommandel(SymmetricTensor{2,3,Tv}, v, offset=offset+NIso+6*i), NKin) # 6 components in SymmetricTensor{2,3}
    return PlasticState(ϵₚ,κ,β)
end
=#

# Note: This could allow symmetric tensors. 
# Could specialize the mandel_stress function to return symmetric tensors
# for isotropic models. 
struct FiniteStrainPlasticResidual{NKin_m1,NIso,T}
    Mred::Tensor{2,3,T,9}
    Δλ::T
    κ::NTuple{NIso, T}
    Mk::NTuple{NKin_m1, Tensor{2,3,T,9}} # Saves Mk[i-1] for i ∈ (2,Nkin), NKin_m1=NKin-1
end

Tensors.get_base(::Type{<:FiniteStrainPlasticResidual{NKin_m1,NIso}}) where{NKin_m1,NIso} = PlasticResidual{NKin_m1,NIso} # needed for frommandel

Tensors.n_components(::Type{<:FiniteStrainPlasticResidual{NKin_m1,NIso}}) where{NKin_m1,NIso} = 10 + NIso + 9*NKin_m1

function Tensors.tomandel!(v::AbstractVector, r::FiniteStrainPlasticResidual{NKin_m1,NIso,T}; offset=0) where {NKin_m1,NIso,T}
    tomandel!(v, r.Mred; offset=offset)
    v[offset+10] = r.Δλ
    v[(offset+11):(offset+10+NIso)] .= r.κ
    for i=1:NKin_m1
        tomandel!(v, r.Mk[i], offset=offset+1+NIso+9*i)
    end
    return v
end

function Tensors.tomandel(r::FiniteStrainPlasticResidual{NKin_m1,NIso,T}) where {NKin_m1,NIso,T}
    v=zeros(T, Tensors.n_components(PlasticResidual{NKin_m1,NIso}))
    return tomandel!(v, r)
end

function Tensors.frommandel(::Type{<:FiniteStrainPlasticResidual{NKin_m1,NIso}}, v::AbstractVector{T}; offset=0) where {T,NKin_m1,NIso}
    Mred = frommandel(Tensor{2,3}, v; offset=offset)
    Δλ = v[offset+10]
    κ = ntuple(i->v[offset+10+i], NIso)
    Mk = ntuple(i->frommandel(Tensor{2,3}, v, offset=offset+1+NIso+9*i), NKin_m1)
    return PlasticResidual(Mred, Δλ, κ, Mk)
end

function MMB.allocate_material_cache(m::FiniteStrainPlastic)
    T = get_base_numbertype(m)
    s = initial_material_state(m)
    F = one(Tensor{2,3})
    x = initial_guess(m, s, F)
    xv = Vector{T}(undef, Tensors.n_components(typeof(x)))
    rf!(r_vector, x_vector) = vector_residual!((x)->residual(x, m, old, F, zero(T)), r_vector, x_vector, x)
    return NewtonCache(xv, rf!)
end

function mandel_stress(m_el::AbstractHyperElastic, Fe::Tensor{2,3})
    Ce = tdot(Fe)
    return 2*Ce ⋅ compute_stress(m_el, Ce)
end

function MMB.material_response(m::Plastic, F::Tensor{2,3}, old::FiniteStrainPlasticState, Δt, cache, extras)
    Fe = F ⋅ inv(old.Fp)
    M = mandel_stress(m.elastic, Fe)
    
    Φ_trial = yield_criterion(m.yield, M - sum(old.Mk), sum(old.κ))

    if Φ_trial < 0
        update_extras!(extras)
        P, dPdF = MMB.material_response(m.elastic, Fe)
        return P, dPdF, old
    else
        x0 = initial_guess(m, old, ϵ)
        rf!(r_vector, x_vector) = vector_residual!(x->residual(x, m, old, ϵ, Δt, cache.resid), r_vector, x_vector, x)
        x_vector = getx(cache.newton)
        tomandel!(x_vector, x0)
        x_vector, ∂R∂X, converged = newtonsolve(x_vector, rf!, cache.newton)
        if converged
            x_sol = frommandel(Tensors.get_base(typeof(x0)), x_vector)
            check_solution(x_sol)
            update_extras!(extras, x_sol, ∂R∂X) # Use dRdx before calling inv!
            # dRdF = 0 = ∂R∂F + ∂R∂X dXdF
            # TODO: Convert to static arrays for speed?
            # TODO: Define get_plastic_result, value_jacobian, and get_plastic_state. 
            ∂R∂F = jacobian(F_ -> tomandel(Vector, residual(x_sol, m, old, frommandel(Tensor{2,3}, F_), Δt)), tomandel(Vector, F))
            dXdF = - ∂R∂X \ ∂R∂F
            P_m, dPdX = value_jacobian(x -> get_plastic_result(m, x, old), x_sol)
            P = frommandel(Tensor{2,3}, P_m)
            dPdF = frommandel(Tensor{4,3}, dPdX * dXdF)
            new = get_plastic_state(m, x_sol, old)
            
            return σ, dσdϵ, new
        else
            throw(MMB.NoLocalConvergence("$(typeof(m)): newtonsolve! didn't converge, ϵ = ", ϵ))
        end
    end
end

check_solution(x::PlasticResidual) = x.Δλ < 0 ? throw(MMB.NoLocalConvergence("Plastic: Invalid solution, x.Δλ = ", x.Δλ, " < 0")) : nothing

# TODO: Could be replaced by exponential map. 
function Fx_time_integration(Fx_old, ν_full, Δλ)
    I2 = one(ν_full)
    return inv(I2 - Δλ * ν_full) ⋅ Fx_old
end

function calculate_neohook_backstress(kh::AbstractKinematicHardening, Fk::Tensor)
    el_law = NeoHooke(get_modulus(kh))
    return mandel_stress(el_law, inv(Fk))
end

# General residual function 
function residual(x::FiniteStrainPlasticResidual, m::FiniteStrainPlastic, old::FiniteStrainPlasticState, F::Tensor, Δt)
    κ = sum(x.κ)
    Φ = yield_criterion(m.yield, x.Mred, κ)

    ν = effective_stress_gradient(m.yield,  x.Mred)
    Fp = Fx_time_integration(old.Fp, ν, x.Δλ)
    Fpinv = inv(Fp)
    Fe = F ⋅ Fpinv
    M = mandel_stress(m.elastic, Fe)

    new_Fk(Fk_old, kh, Mk) = Fx_time_integration(Fk_old, get_evolution(kh, ν, Mk), x.Δλ)

    Mk1 = M - (x.Mred + sum(x.Mk))
    Mk_guess = tuple(Mk1, x.Mk...)
    Fk = map(new_Fk, old.Fk, m.kinematic, Mk_guess)

    Mk = map(calculate_neohook_backstress, m.kinematic, Fk)
    
    RM = x.Mred - (M - sum(Mk))
    Rλ = yield_residual(m.overstress, Φ, x.Δλ, Δt, initial_yield_limit(m.yield) + κ)
    Rκ = map((ih, κ, κold) -> κ - κold - x.Δλ * get_evolution(ih, κ), m.isotropic, x.κ, old.κ)
    RMk = map(-, Mk[2:end], x.Mk)
    
    return FiniteStrainPlasticResidual(RM, Rλ, Rκ, RMk)
end

function initial_guess(m::FiniteStrainPlastic, old::FiniteStrainPlasticState, M_trial)
    Mk = map(calculate_neohook_backstress, m.kinematic, old.Fk)
    Δλ = 0.0
    return FiniteStrainPlasticResidual(M_trial-sum(Mk), Δλ, old.κ, Mk[2:end])
end

#=
# Functions for conversion between material and parameter vectors
function MMB.get_num_params(m::Plastic)
    return sum((
        get_num_params(m.elastic),
        get_num_params(m.yield),
        sum(get_num_params, m.isotropic),
        sum(get_num_params, m.kinematic),
        get_num_params(m.overstress)))
end

function vector2materialtuple(v::AbstractVector, materialtuple; offset=0)
    n = 0
    mout = map(materialtuple) do mt
        m = vector2material(v, mt; offset=(offset+n))
        n += get_num_params(m)
        m
    end
    return mout, n
end

function MMB.vector2material(v::AbstractVector, m::Plastic; offset=0)
    i = offset
    elastic = vector2material(v, m.elastic, offset=i); i += get_num_params(m.elastic)
    yield = vector2material(v, m.yield, offset=i); i += get_num_params(m.yield)
    
    isotropic, nisoparam = vector2materialtuple(v, m.isotropic, offset=i); i+=nisoparam
    kinematic, nkinparam = vector2materialtuple(v, m.kinematic, offset=i); i+=nkinparam

    overstress = vector2material(v, m.overstress, offset=i); # i += get_num_params(m.overstress)

    return Plastic(;elastic, yield, isotropic, kinematic, overstress)
end

function MMB.material2vector!(v::AbstractVector, m::Plastic; offset=0)
    i = offset
    material2vector!(v, m.elastic; offset=i); i += get_num_params(m.elastic)
    material2vector!(v, m.yield; offset=i); i += get_num_params(m.yield)
    for iso in m.isotropic
        material2vector!(v, iso;offset=i); i+= get_num_params(iso)
    end
    for kin in m.kinematic
        material2vector!(v, kin; offset=i); i+= get_num_params(kin)
    end
    material2vector!(v, m.overstress; offset=i); # i += get_num_params(m.overstress)
    return v
end

# Show method
function Base.show(io::IO, M::MIME"text/plain", m::Plastic)
    println(io, "Plastic material with")
    show(io, M, m.elastic)
    show(io, M, m.yield)
    if length(m.kinematic) == 1
        print(io, "Kinematic hardening: ")
        show(io, M, m.kinematic[1])
    else
        println(io, "Kinematic hardening laws:")
        foreach(kin->show(io, M, kin), m.kinematic)
    end
    if length(m.kinematic) == 1
        print(io, "Isotropic hardening: ")
        show(io, M, m.isotropic[1])
    else    
        println(io, "Isotropic hardening laws:")
        foreach(iso->show(io, M, iso), m.isotropic)
    end
    show(io, M, m.overstress)
end
=#