"""
Plastic(;elastic, yield, isotropic, kinematic, overstress)

A plasticity model with modular elastic law, yield criterion, 
multiple modular isotropic and kinematic hardening contributions,
and either rate-independent or viscoplastic response.
# Arguments
- `elastic::AbstractMaterial`: Elastic law, see e.g. [`LinearElastic`](@ref)
- `yield::YieldCriterion`: Yield criterion, including the initial yield limit. If yield::Real is given, `VonMises(yield)` is used. 
- `isotropic::Union{AbstractIsotropicHardening,Tuple}`: Isotropic hardening laws, see e.g. [`Voce`](@ref)
- `kinematic::Union{AbstractKinematicHardening,Tuple}`: Kinematic hardening laws, see e.g. [`ArmstrongFrederick`](@ref)
- `overstress`::Union{RateIndependent,OverstressFunction}: Rate dependence, see e.g. [`NortonOverstress`](@ref)

The model response is given by the laws in `elastic`, `yield`, `isotropic`, `kinematic`, and `overstress`.
The generic model equations are described below. 

The stress is calculated from the elastic strains, ``\\boldsymbol{\\epsilon}_\\mathrm{e}``, obtained via the 
additive decomposition, ``\\boldsymbol{\\epsilon} = \\boldsymbol{\\epsilon}_\\mathrm{e} + \\boldsymbol{\\epsilon}_\\mathrm{p}``. 
The elastic law is specified by `m.elastic` and is evaluated by giving it the elastic strain. 

A yield criterion of the type 
```math
\\varPhi = f\\left( \\boldsymbol{\\sigma} - \\boldsymbol{\\beta} \\right) - \\left[Y_0 - \\kappa\\right]
```
is assumed. Here, ``\\boldsymbol{\\beta} = \\sum_{i=1}^{N_\\mathrm{kin}} \\boldsymbol{\\beta}_i`` is the total back-stress, 
and ``\\kappa = \\sum_{i=1}^{N_\\mathrm{iso}} \\kappa_i`` is the total isotropic hardening stress. The initial yield limit 
is passed to the yield criterion along with potentially other parameters. 
The evolution laws for ``\\boldsymbol{\\beta}_i`` and ``\\kappa_i`` are given by the kinematic and isotropic hardening laws.

Associative plastic flow is used to obtain the plastic strains,
```math
\\dot{\\epsilon}_{\\mathrm{p}} = \\dot{\\lambda} \\left.\\frac{\\partial f}{\\partial \\boldsymbol{\\sigma}}\\right\\vert_{\\left( \\boldsymbol{\\sigma} - \\boldsymbol{\\beta} \\right)}
= \\dot{\\lambda} \\boldsymbol{\\nu}
```

The isotropic hardening is formulated as
```math
\\kappa = \\sum_{i=1}^{N_{\\mathrm{iso}}} g_{\\mathrm{iso},i}(\\lambda)
```
where ``g_{\\mathrm{iso},i}(\\lambda)`` is specified by `m.isotropic[i]` (see [Isotropic hardening](@ref))

Kinematic hardening is formulated as
```math
\\dot{\\boldsymbol{\\beta}}_i = \\dot{\\lambda} g_{\\mathrm{kin},i}(\\nu, \\boldsymbol{\\beta}_i)
```
where ``g_{\\mathrm{kin},i}(\\boldsymbol{\\nu}, \\boldsymbol{\\beta}_i)`` is specified by `m.kinematic[i]`
and ``i\\in[1,N_\\mathrm{kin}]`` (see [Kinematic hardening](@ref))]

If `overstress=RateIndependent()`, the plastic multiplier, ``\\lambda``, is obtained via the KKT-conditions,
```math
\\dot{\\lambda} \\geq 0, \\quad \\varPhi \\leq 0, \\quad \\dot{\\lambda}\\varPhi = 0
```
Otherwise, the overstress function, ``\\eta(\\varPhi)``, determines the evolution of ``\\lambda`` as 
```math
\\dot{\\lambda} = \\eta(\\varPhi, (Y_0 + \\kappa))
```

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

get_base_numbertype(m::Plastic) = get_base_numbertype(m.elastic)

# Definition of material state
struct PlasticState{NKin,NIso,Tϵₚ,Tκ,Tβ} <: AbstractMaterialState
    # Useful to have different types here?
    ϵₚ::SymmetricTensor{2,3,Tϵₚ,6}
    κ::NTuple{NIso, Tκ}
    β::NTuple{NKin, SymmetricTensor{2,3,Tβ,6}}
end

function MMB.initial_material_state(m::Plastic)
    T = get_base_numbertype(m)
    PlasticState(zero(SymmetricTensor{2,3,T}), 
                 ntuple(i->get_initial_value(m.isotropic[i]), length(m.isotropic)), 
                 ntuple(i->zero(SymmetricTensor{2,3,T}), length(m.kinematic)))
end

Tensors.n_components(::Type{<:PlasticState{NK, NI}}) where {NK, NI} = NI + 6*(1+NK)

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
    T = get_base_numbertype(m)
    s = initial_material_state(m)
    ϵ = zero(SymmetricTensor{2,3})
    x = initial_guess(m, s, ϵ)
    xv = Vector{T}(undef, Tensors.n_components(typeof(x)))
    rf!(r_vector, x_vector) = vector_residual!((x)->residual(x, m, old, ϵ, zero(T), residual_cache), r_vector, x_vector, x)
    return NewtonCache(xv, rf!)
end

function MMB.allocate_material_cache(m::Plastic)
    return PlasticCache(get_newton_cache(m, nothing), nothing)
end

function elastic_response(m::Plastic, ϵ::SymmetricTensor{2,3}, old::AbstractMaterialState, Δt=nothing)
    σ, dσdϵ, _ = material_response(m.elastic, ϵ-old.ϵₚ, initial_material_state(m.elastic))
    return σ,dσdϵ
end

function iso_yield_limit(m::Plastic, state::Union{AbstractMaterialState, PlasticResidual})
    minyield = 40.0 # Bound hardcoded here for convenience, should be lifted out for better flexibility later
    return smoothbound(m.σ_y0 + sum(state.κ), minyield, minyield/2)
end

function MMB.material_response(m::Plastic, ϵ::SymmetricTensor{2,3}, old::PlasticState{NKin,NIso}, Δt=nothing, cache=allocate_material_cache(m), extras::AbstractExtraOutput=NoExtraOutput(); options = nothing) where {NKin,NIso}

    σ_trial, dσdϵ_elastic = elastic_response(m, ϵ, old)
    Φ_trial = yield_criterion(m.yield, σ_trial-sum(old.β), sum(old.κ))

    if Φ_trial < 0
        update_extras!(extras)
        return σ_trial, dσdϵ_elastic, old
    else
        x = initial_guess(m, old, ϵ)
        rf!(r_vector, x_vector) = vector_residual!((x)->residual(x, m, old, ϵ, Δt, cache.resid), r_vector, x_vector, x)
        x_vector = getx(cache.newton); tomandel!(x_vector, x)
        x_vector, dRdx, converged = newtonsolve(x_vector, rf!, cache.newton)
        if converged
            x = frommandel(PlasticResidual{NKin,NIso}, x_vector)
            check_solution(x)
            inv_J_σσ = frommandel(SymmetricTensor{4,3}, inv(dRdx))
            typeof(m.elastic) <: LinearElastic || error("Only LinearElastic elasticity supported") # Otherwise, the following expression is not true
            dσdϵ = inv_J_σσ ⊡ dσdϵ_elastic
            σ, new = get_plastic_result(x, old) # In paper denote this as f_y(X(ϵ,ⁿs,p), ϵ, ⁿs, p)
                                                # But in this implementation, only f_y(X(ϵ,ⁿs,p), ⁿs)
            update_extras!(extras, x, dRdx)
            return σ, dσdϵ, new
        else
            @show ϵ
            throw(MMB.NoLocalConvergence("$(typeof(m)): newtonsolve! didn't converge"))
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
    Rκ = ntuple(i-> x.κ[i] - old.κ[i] - x.Δλ * get_evolution(m.isotropic[i], x.κ[i]), NIso)
    Rβ = ntuple((i) -> x.β[i] - old.β[i] - x.Δλ * get_evolution(m.kinematic[i], ν, x.β[i]), NKin)

    return PlasticResidual(Rσ, Rλ, Rκ, Rβ)
end

function initial_guess(m::Plastic, old::PlasticState, ϵ)
    σ_trial = calculate_stress(m.elastic, ϵ-old.ϵₚ)
    Δλ = 0.0
    return PlasticResidual(σ_trial, Δλ, old.κ, old.β)
end

function get_plastic_result(x::PlasticResidual, old::PlasticState)
    σ_red_dev = dev(x.σ) - sum(x.β)
    ϵₚ = calculate_plastic_strain(old, σ_red_dev * ((3/2)/vonmises(σ_red_dev)), x.Δλ)
    return x.σ, PlasticState(ϵₚ, x.κ, x.β)
end

function calculate_elastic_strain(old::PlasticState, ϵ, ν, Δλ)
    return ϵ - calculate_plastic_strain(old, ν, Δλ)
end

function calculate_plastic_strain(old::PlasticState, ν, Δλ)
    return old.ϵₚ + Δλ*ν
end

effective_vm(σ::AbstractTensor, β::Tuple) = effective_vm(σ, sum(β))
function effective_vm(σ::AbstractTensor, β::AbstractTensor)
    σ_red_dev = dev(σ) - β
    σ_vm = vonmises(σ_red_dev)
    ν = σ_red_dev * ((3/2)/σ_vm) 
    return σ_vm, ν
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
