"""
    Plastic(;elastic, yield, isotropic, kinematic, overstress)

A small-strain plasticity model with modular elastic laws, yield criteria, 
multiple isotropic and kinematic hardening contributions,
and either rate-independent or viscoplastic response.
# Keyword arguments
- `elastic::AbstractMaterial`\\
  Elastic law, see e.g. [`LinearElastic`](@ref)
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
m = Plastic(elastic = LinearElastic(E=210.e3, ν=0.3),
            yield = 100.0,
            isotropic = (Voce(Hiso=-100.e3, κ∞=-100.0),Voce(Hiso=10.e3, κ∞=200.0)),
            kinematic = (ArmstrongFrederick(Hkin=200.e3, β∞=200.0),
                         OhnoWang(Hkin=1000.e3, β∞=200.0, m=3.0)),
            overstress = NortonOverstress(;tstar=1.0, nexp=2.0))
```

"""
struct Plastic{E,YLD,IsoH,KinH,OS} <: AbstractMaterial
    elastic::E          # Elastic definition
    yield::YLD          # Initial yield limit
    isotropic::IsoH     # Tuple of isotropic hardening definitions
    kinematic::KinH     # Tuple of kinematic hardening definitions
    overstress::OS      # Overstress function
end
function Plastic(;elastic, yield, isotropic=nothing, kinematic=nothing, overstress=RateIndependent())
    iso = maketuple_or_nothing(isotropic)
    kin = maketuple_or_nothing(kinematic)
    yld = default_yield_criteria(yield)
    return Plastic(elastic, yld, iso, kin, overstress)
end

MMB.get_parameter_type(m::Plastic) = MMB.get_parameter_type(m.elastic)

# Definition of material state
struct PlasticState{NKin,NIso,Tϵp,Tκ,Tβ} <: AbstractMaterialState
    # Useful to have different types here?
    ϵp::SymmetricTensor{2,3,Tϵp,6}
    κ::NTuple{NIso, Tκ}
    β::NTuple{NKin, SymmetricTensor{2,3,Tβ,6}}
end

function MMB.initial_material_state(m::Plastic)
    T = MMB.get_parameter_type(m)
    PlasticState(zero(SymmetricTensor{2,3,T}), 
                 ntuple(i->get_initial_value(m.isotropic[i]), length(m.isotropic)), 
                 ntuple(i->zero(SymmetricTensor{2,3,T}), length(m.kinematic)))
end

Tensors.n_components(::Type{<:PlasticState{NK, NI}}) where {NK, NI} = NI + 6*(1+NK)

MMB.get_num_statevars(m::Plastic) = length(m.isotropic) + 6*(1 + length(m.kinematic))

function Tensors.tomandel!(v::AbstractVector, s::PlasticState{NKin,NIso,Tϵp,Tκ,Tβ}; offset=0) where {NKin,NIso,Tϵp,Tκ,Tβ}
    tomandel!(v, s.ϵp; offset=offset)
    v[(offset+7):(offset+6+NIso)] .= s.κ
    for i=1:NKin
        tomandel!(v, s.β[i], offset=offset+NIso+6*i) # 6 components in SymmetricTensor{2,3}
    end
    return v
end
function Tensors.tomandel(s::PlasticState{NKin,NIso,Tϵp,Tκ,Tβ}) where {NKin,NIso,Tϵp,Tκ,Tβ}
    v = zeros(promote_type(Tϵp,Tκ,Tβ), Tensors.n_components(PlasticState{NKin,NIso}))
    return tomandel!(v, s)
end

function Tensors.frommandel(::Type{<:PlasticState{NKin,NIso}}, v::AbstractVector{Tv}; offset=0) where {Tv,NKin,NIso}
    ϵp = frommandel(SymmetricTensor{2,3}, v; offset=offset)
    κ = ntuple(i->v[offset+6+i], NIso)
    β = ntuple(i->frommandel(SymmetricTensor{2,3,Tv}, v, offset=offset+NIso+6*i), NKin) # 6 components in SymmetricTensor{2,3}
    return PlasticState(ϵp,κ,β)
end

struct PlasticResidual{NKin,NIso,Tσ,Tλ,Tκ,Tβ}
    σ::SymmetricTensor{2,3,Tσ,6}
    Δλ::Tλ
    κ::NTuple{NIso, Tκ}
    β::NTuple{NKin, SymmetricTensor{2,3,Tβ,6}}
end

Tensors.get_base(::Type{<:PlasticResidual{NKin,NIso}}) where{NKin,NIso} = PlasticResidual{NKin,NIso} # needed for frommandel

Tensors.n_components(::Type{<:PlasticResidual{NKin,NIso}}) where{NKin,NIso} = 7 + NIso + 6*NKin

function Tensors.tomandel!(v::AbstractVector, r::PlasticResidual{NKin,NIso,Tσ,Tλ,Tκ,Tβ}; offset=0) where {NKin,NIso,Tσ,Tλ,Tκ,Tβ}
    tomandel!(v, r.σ; offset=offset)
    v[offset+7] = r.Δλ
    v[(offset+8):(offset+7+NIso)] .= r.κ
    for i=1:NKin
        tomandel!(v, r.β[i], offset=offset+1+NIso+6*i) # 6 components in SymmetricTensor{2,3}
    end
    return v
end

function Tensors.tomandel(r::PlasticResidual{NKin,NIso,Tσ,Tλ,Tκ,Tβ}) where {NKin,NIso,Tσ,Tλ,Tκ,Tβ}
    v=zeros(promote_type(Tσ,Tλ,Tκ,Tβ), Tensors.n_components(PlasticResidual{NKin,NIso}))
    return tomandel!(v, r)
end

function Tensors.frommandel(::Type{<:PlasticResidual{NKin,NIso}}, v::AbstractVector{Tv}; offset=0) where {Tv,NKin,NIso}
    σ = frommandel(SymmetricTensor{2,3}, v; offset=offset)
    Δλ = v[offset+7]
    κ = ntuple(i->v[offset+7+i], NIso)
    β = ntuple(i->frommandel(SymmetricTensor{2,3,Tv}, v, offset=offset+1+NIso+6*i), NKin)
    return PlasticResidual(σ,Δλ,κ,β)
end

struct PlasticCache{TN,TR}
    newton::TN  # NewtonCache if needed, otherwise nothing
    resid::TR   # Cache used in residual function if needed, otherwise nothing
end

function get_newton_cache(m::Plastic, residual_cache)
    T = MMB.get_parameter_type(m)
    s = initial_material_state(m)
    ϵ = zero(SymmetricTensor{2,3})
    x = initial_guess(m, s, ϵ)
    xv = Vector{T}(undef, Tensors.n_components(typeof(x)))
    rf!(r_vector, x_vector) = vector_residual!((x)->residual(x, m, old, ϵ, zero(T), residual_cache), r_vector, x_vector, x)
    return NewtonCache(xv)
end

function MMB.allocate_material_cache(m::Plastic)
    return PlasticCache(get_newton_cache(m, nothing), nothing)
end

function elastic_response(m::Plastic, ϵ::SymmetricTensor{2,3}, old::AbstractMaterialState, Δt=nothing)
    σ, dσdϵ, _ = material_response(m.elastic, ϵ-old.ϵp, initial_material_state(m.elastic))
    return σ,dσdϵ
end

function MMB.material_response(m::Plastic, ϵ::SymmetricTensor{2,3}, old::PlasticState{NKin,NIso}, Δt, cache, extras) where {NKin,NIso}

    σ_trial, dσdϵ_elastic = elastic_response(m, ϵ, old)
    Φ_trial = yield_criterion(m.yield, σ_trial-sum(old.β), sum(old.κ))

    if Φ_trial < 0
        update_extras!(extras)
        return σ_trial, dσdϵ_elastic, old
    else
        x = initial_guess(m, old, ϵ)
        rf!(r_vector, x_vector) = vector_residual!(xx->residual(xx, m, old, ϵ, Δt, cache.resid), r_vector, x_vector, x)
        x_vector = getx(cache.newton)
        tomandel!(x_vector, x)
        x_vector, dRdx, converged = newtonsolve(rf!, x_vector, cache.newton)
        if converged
            x_sol = frommandel(PlasticResidual{NKin,NIso}, x_vector)
            check_solution(x_sol)
            inv_dRdx = Newton.inv!(dRdx, cache.newton)
            inv_J_σσ = frommandel(SymmetricTensor{4,3}, inv_dRdx)
            typeof(m.elastic) <: LinearElastic || error("Only LinearElastic elasticity supported") # Otherwise, the following expression is not true
            dσdϵ = inv_J_σσ ⊡ dσdϵ_elastic
            σ, new = get_plastic_result(m, x_sol, old) # In general, f_y(X(ϵ,ⁿs,p), ϵ, ⁿs, p)
                                                       # But here,   f_y(X(ϵ,ⁿs,p), ⁿs, p), suffices
            update_extras!(extras, x_sol, inv_dRdx)
            return σ, dσdϵ, new
        else
            throw(MMB.NoLocalConvergence("$(typeof(m)): newtonsolve! didn't converge, ϵ = ", ϵ))
        end
    end
end

check_solution(x::PlasticResidual) = x.Δλ < 0 ? throw(MMB.NoLocalConvergence("Plastic: Invalid solution, x.Δλ = ", x.Δλ, " < 0")) : nothing

# General residual function 
function residual(x::PlasticResidual{NKin,NIso}, m::Plastic, old::PlasticState, ϵ, Δt, cache) where{NKin,NIso}
    σ_red = x.σ - sum(x.β)
    ν = effective_stress_gradient(m.yield,  σ_red)
    ϵₑ = calculate_elastic_strain(old, ϵ, ν, x.Δλ)    # Using assumption of associative plasticity

    κ = sum(x.κ)
    Φ = yield_criterion(m.yield, σ_red, κ)
    σy = initial_yield_limit(m.yield) + κ

    Rσ = x.σ - calculate_stress(m.elastic, ϵₑ)       # Using the specific elastic law
    Rλ = yield_residual(m.overstress, Φ, x.Δλ, Δt, σy)
    Rκ = map((ih, κ, κold) -> κ - κold - x.Δλ * get_evolution(ih, κ), m.isotropic, x.κ, old.κ)
    Rβ_fun(kh, β, βold) = β - βold + (x.Δλ * 2*get_modulus(kh)/3) * get_evolution(kh, ν, β)
    Rβ = map(Rβ_fun, m.kinematic, x.β, old.β)
    
    return PlasticResidual(Rσ, Rλ, Rκ, Rβ)
end

function initial_guess(m::Plastic, old::PlasticState, ϵ)
    σ_trial = calculate_stress(m.elastic, ϵ-old.ϵp)
    Δλ = 0.0
    return PlasticResidual(σ_trial, Δλ, old.κ, old.β)
end

function get_plastic_result(m::Plastic, x::PlasticResidual, old::PlasticState)
    ν = effective_stress_gradient(m.yield,  x.σ - sum(x.β))
    ϵp = old.ϵp + x.Δλ*ν
    return x.σ, PlasticState(ϵp, x.κ, x.β)
end

function calculate_elastic_strain(old::PlasticState, ϵ, ν, Δλ)
    return ϵ - (old.ϵp + Δλ*ν)
end

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
    n = Ref(0)
    mout = map(materialtuple) do mt
        m = vector2material(v, mt; offset=(offset+n[]))
        n[] += get_num_params(mt)
        m
    end
    return mout, n[]
end

function MMB.vector2material(v::AbstractVector, m::Plastic; offset=0)
    i = offset
    elastic = vector2material(v, m.elastic, offset=i); i += get_num_params(m.elastic)
    yield = vector2material(v, m.yield, offset=i); i += get_num_params(m.yield)
    
    isotropic, nisoparam = vector2materialtuple(v, m.isotropic, offset=i)
    i += nisoparam
    kinematic, nkinparam = vector2materialtuple(v, m.kinematic, offset=i)
    i += nkinparam

    overstress = vector2material(v, m.overstress, offset=i); # i += get_num_params(m.overstress)
    # The following constructor call doesn't seem to be type stable when used in e.g. differentiate_material. 
    # Should be checked for (a) v and m same eltype and (b) v Dual and m Float64 eltypes.
    #return Plastic(;elastic, yield, isotropic, kinematic, overstress)
    return Plastic(elastic, yield, isotropic, kinematic, overstress)
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
    if length(m.isotropic) == 1
        print(io, "Isotropic hardening: ")
        show(io, M, m.isotropic[1])
    else    
        println(io, "Isotropic hardening laws:")
        foreach(iso->show(io, M, iso), m.isotropic)
    end
    show(io, M, m.overstress)
end
