"""
    FiniteStrainPlastic(;elastic, yield, isotropic, kinematic, overstress)

A finite-strain plasticity model with modular elastic laws, yield criteria, 
multiple isotropic and kinematic hardening contributions,
and either rate-independent or viscoplastic response.
# Keyword arguments
- `elastic::AbstractHyperElastic`\\
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
m = FiniteStrainPlastic(elastic = CompressibleNeoHooke(G=80.e3, K=160.e3),
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


MMB.get_parameter_type(m::FiniteStrainPlastic) = typeof(initial_yield_limit(m.yield))

# Definition of material state
struct FiniteStrainPlasticState{NKin,NIso,TFp,Tκ<:NTuple{NIso},TFk<:NTuple{NKin}} <: AbstractMaterialState
    Fp::TFp
    κ::Tκ
    Fk::TFk
end

function MMB.initial_material_state(m::FiniteStrainPlastic)
    T = MMB.get_parameter_type(m)
    I2 = one(Tensor{2,3,T})
    FiniteStrainPlasticState(
        I2,                                  # Fp
        map(get_initial_value, m.isotropic), # κ
        map(Returns(I2), m.kinematic))       # Fk
end

# Note: This could allow symmetric tensors. 
# Could specialize the mandel_stress function to return symmetric tensors
# for isotropic models. 
struct FiniteStrainPlasticResidual{NKin_m1,NIso,TM,Tλ,Tκ<:NTuple{NIso},TMk<:NTuple{NKin_m1}}
    Mred::TM
    Δλ::Tλ
    κ::Tκ
    Mk::TMk # Saves Mk[i-1] for i ∈ (2,Nkin), NKin_m1=NKin-1
end
function FiniteStrainPlasticResidual(Mred, Δλ, κ::Tuple, Mk::Tuple)
    Tκ = get_promoted_type(κ...)
    TMk = get_promoted_type(Mk...)
    return FiniteStrainPlasticResidual(Mred, Δλ, map(x->convert(Tκ, x), κ), map(x->convert(TMk, x), Mk))
end

Tensors.get_base(::Type{<:FiniteStrainPlasticResidual{NKin_m1,NIso}}) where{NKin_m1,NIso} = FiniteStrainPlasticResidual{NKin_m1,NIso} # needed for frommandel

Tensors.n_components(::Type{<:FiniteStrainPlasticResidual{NKin_m1,NIso}}) where{NKin_m1,NIso} = 10 + NIso + 9*NKin_m1

function Tensors.tomandel(::Type{<:SVector}, r::FiniteStrainPlasticResidual{NKin_m1,NIso}) where {NKin_m1,NIso}
    Ms = tomandel(SVector, r.Mred)
    Δλ = r.Δλ
    κs = SVector(r.κ)
    Mks = map(x->tomandel(SVector, x), r.Mk)
    return static_vector(Ms, Δλ, κs, Mks...)
end

function Tensors.frommandel(::Type{<:FiniteStrainPlasticResidual{NKin_m1,NIso}}, v::AbstractVector{T}; offset=0) where {T,NKin_m1,NIso}
    Mred = frommandel(Tensor{2,3}, v; offset=offset)
    Δλ = v[offset+10]
    κ = ntuple(i->v[offset+10+i], NIso)
    Mk = ntuple(i->frommandel(Tensor{2,3}, v, offset=offset+1+NIso+9*i), NKin_m1)
    return FiniteStrainPlasticResidual(Mred, Δλ, κ, Mk)
end

function MMB.allocate_material_cache(m::FiniteStrainPlastic)
    T = MMB.get_parameter_type(m)
    s = initial_material_state(m)
    F = one(Tensor{2,3})
    x = initial_guess(m, s, F)
    xv = Vector{T}(undef, Tensors.n_components(typeof(x)))
    rf!(r_vector, x_vector) = vector_residual!((x)->residual(x, m, old, F, zero(T)), r_vector, x_vector, x)
    return NewtonCache(xv)
end

function mandel_stress(m_el::AbstractHyperElastic, Fe::Tensor{2,3})
    Ce = tdot(Fe)
    return Ce ⋅ compute_stress(m_el, Ce)
end

function MMB.material_response(m::FiniteStrainPlastic, F::Tensor{2,3}, old::FiniteStrainPlasticState, Δt, cache, extras)
    Fpinv = inv(old.Fp)
    Fe = F ⋅ Fpinv
    M = mandel_stress(m.elastic, Fe)
    Mk_old = map(calculate_neohook_backstress, m.kinematic, old.Fk)
    Φ_trial = yield_criterion(m.yield, - sum(Mk_old; init=-M), sum(old.κ; init=0))

    if Φ_trial < 0
        update_extras!(extras)
        dPdF, P = gradient(F_->calculate_PKstress(m, old, F_), F, :all)
        return P, dPdF, old
    else
        x0 = initial_guess(m, old, M)
        rf(x) = tomandel(SVector, residual(frommandel(typeof(x0), x), m, old, F, Δt))
        x_vector, ∂R∂X, converged = newtonsolve(rf, tomandel(SVector, x0))

        if converged
            x_sol = frommandel(typeof(x0), x_vector)
            check_solution(x_sol)
            update_extras!(extras, x_sol, ∂R∂X) # Use dRdx before calling inv!
            # dPdF = ∂P∂F + ∂P∂X dXdF
            # dRdF = 0 = ∂R∂F + ∂R∂X dXdF => dXdF = - ∂R∂X\∂R∂F
            
            F_R_fun = Tensor2VectorFun(F, F_-> residual(x_sol, m, old, F_, Δt))
            ∂R∂F = ForwardDiff.jacobian(F_R_fun, tomandel(SVector, F))
            
            dXdF = - ∂R∂X \ ∂R∂F

            P = calculate_PKstress(m, x_sol, old, F)
            P_X_fun = Tensor2VectorFun(x_sol, x -> calculate_PKstress(m, x, old, F))
            ∂P∂X = ForwardDiff.jacobian(P_X_fun, tomandel(SVector, x_sol))

            ∂P∂F = gradient(F_->calculate_PKstress(m, x_sol, old, F_), F)

            dPdF = ∂P∂F + frommandel(Tensor{4,3}, ∂P∂X * dXdF)
            new, _ = get_plastic_state(x_sol, m, old, F, Δt)
            
            return P, dPdF, new
        else
            throw(MMB.NoLocalConvergence("$(typeof(m)): newtonsolve! didn't converge, F = ", F))
        end
    end
end

struct Tensor2VectorFun{T, F<:Function}
    tt::Type{T}
    f::F
end
Tensor2VectorFun(::T, f::Function) where T = Tensor2VectorFun(Tensors.get_base(T), f)
#(tvf::Tensor2VectorFun{T})(x::Vector) where T = tomandel(tvf.f(frommandel(T, x)))
(tvf::Tensor2VectorFun{T})(x::SVector) where T = tomandel(SVector, tvf.f(frommandel(T, x)))

function get_plastic_state(x::FiniteStrainPlasticResidual, m::FiniteStrainPlastic, old::FiniteStrainPlasticState, F::Tensor, Δt)
    ν = effective_stress_gradient(m.yield,  x.Mred)
    Fp = Fx_time_integration(old.Fp, ν, x.Δλ)
    Fpinv = inv(Fp)
    Fe = F ⋅ Fpinv
    M = mandel_stress(m.elastic, Fe)

    Mk1 = M - sum(x.Mk; init=x.Mred)
    T_Mk = promote_type(typeof(Mk1), eltype(x.Mk))
    
    new_Fk(Fk_old, kh, Mk) = Fx_time_integration(Fk_old, get_evolution(kh, ν, Mk), x.Δλ)
    Mk_guess = tuple(Mk1, x.Mk...)
    Fk = map(new_Fk, old.Fk, m.kinematic, Mk_guess)
    new = FiniteStrainPlasticState(Fp, x.κ, Fk)    

    return new, M
end

function calculate_PKstress(m::FiniteStrainPlastic, state::FiniteStrainPlasticState, F::Tensor)
    return calculate_PKstress(m, state.Fp, F)
end
function calculate_PKstress(m::FiniteStrainPlastic, x::FiniteStrainPlasticResidual, old::FiniteStrainPlasticState, F::Tensor)
    ν = effective_stress_gradient(m.yield,  x.Mred)
    Fp = Fx_time_integration(old.Fp, ν, x.Δλ)
    return calculate_PKstress(m, Fp, F)
end
function calculate_PKstress(m::FiniteStrainPlastic, Fp::Tensor, F::Tensor)
    Fpinv = inv(Fp)
    Fe = F ⋅ Fpinv
    Ce = tdot(Fe)
    Se = compute_stress(m.elastic, Ce)
    P = Fe ⋅ Se ⋅ transpose(Fpinv)
    return P
end

check_solution(x::FiniteStrainPlasticResidual) = x.Δλ < 0 ? throw(MMB.NoLocalConvergence("Plastic: Invalid solution, x.Δλ = ", x.Δλ, " < 0")) : nothing

# TODO: Could be replaced by exponential map. 
function Fx_time_integration(Fx_old, ν_full, Δλ)
    I2 = one(ν_full)
    return inv(I2 - Δλ * ν_full) ⋅ Fx_old
end

function calculate_neohook_backstress(kh::AbstractKinematicHardening, Fk::Tensor)
    el_law = NeoHooke(get_modulus(kh)/3)
    return mandel_stress(el_law, inv(Fk))
end

# General residual function 
function residual(x::FiniteStrainPlasticResidual, m::FiniteStrainPlastic, old::FiniteStrainPlasticState, F::Tensor, Δt)
    κ = sum(x.κ; init=0)
    new, M = get_plastic_state(x, m, old, F, Δt)
    Fk = new.Fk
    Φ = yield_criterion(m.yield, x.Mred, κ)

    Mk = map(calculate_neohook_backstress, m.kinematic, Fk)
    
    RM = x.Mred + sum(Mk; init=-M)
    Rλ = yield_residual(m.overstress, Φ, x.Δλ, Δt, initial_yield_limit(m.yield) + κ)
    Rκ = map((ih, κ, κold) -> κ - κold - x.Δλ * get_evolution(ih, κ), m.isotropic, x.κ, old.κ)
    RMk = map(-, Mk[2:end], x.Mk)
    
    return FiniteStrainPlasticResidual(RM, Rλ, Rκ, RMk)
end

function initial_guess(m::FiniteStrainPlastic, old::FiniteStrainPlasticState, M_trial)
    Mk = map(calculate_neohook_backstress, m.kinematic, old.Fk)
    Δλ = 0.0
    return FiniteStrainPlasticResidual(-sum(Mk; init=-M_trial), Δλ, old.κ, Mk[2:end])
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