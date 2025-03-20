"""
    Maxwell(; G, η, t)

Create a maxwell chain element with shear stiffness `G`. Either the viscosity, `η`, 
or the relaxation time, `t = η / G`, should be supplied (not both).
"""
struct Maxwell{T}
    G::T
    η::T
end
function Maxwell(; G, η = nothing, t = nothing)
    if η === nothing
        t === nothing && throw(ArgumentError("Either t or η must be given, but not both"))
        return Maxwell(G, G * t)
    else
        t !== nothing && throw(ArgumentError("Either t or η must be given, but not both"))
        return Maxwell(G, η)
    end
end

function calculate_viscous_strain(m::Maxwell, ϵ::SymmetricTensor, ϵv_old::SymmetricTensor, Δt::Number)
    # ϵv - ϵv_old = Δt * σv / m.η
    # σv = 2G * (ϵdev - ϵv)
    # ϵv * (1 + Δt * 2G / m.η) = ϵv_old +  Δt * 2G * ϵdev / m.η
    Δt_2G_div_η = Δt * 2 * m.G / m.η
    return (ϵv_old +  Δt_2G_div_η * dev(ϵ)) / (1 + Δt_2G_div_η)
end

function calculate_stress(m::Maxwell, ϵ::SymmetricTensor, ϵv_old::SymmetricTensor, Δt::Number)
    ϵv = calculate_viscous_strain(m, ϵ, ϵv_old, Δt)
    return 2 * m.G * (dev(ϵ) - ϵv)
end

"""
    GeneralizedMaxwell(elastic::LinearElastic, chains::Maxwell...)

Create a generalized Maxwell model with an arbitrary number of Maxwell chains.
The `elastic` part refers to the long-term stiffness contribution. 

!!! note "WIP"
    The `GeneralizedMaxwell` and `Maxwell` are currently quite specific to isotropic behavior,
    this may change in the future with related changes to the constructors.

"""
struct GeneralizedMaxwell{ET, T, num} <: AbstractMaterial
    base::ET
    chains::NTuple{num, Maxwell{T}}
end

function GeneralizedMaxwell(base::LinearElastic, chains::Maxwell...)
    return GeneralizedMaxwell(base, chains)
end

struct GeneralizedMaxwellState{T, num} <: AbstractMaterialState
    ϵv::NTuple{num, SymmetricTensor{2,3,T,6}}
end

function MMB.initial_material_state(::GeneralizedMaxwell{<:Any, <:Any, num}) where {num}
    T = Float64 # temp
    return GeneralizedMaxwellState(ntuple(_ -> zero(SymmetricTensor{2, 3, T}), num))
end

function calculate_stress(m::GeneralizedMaxwell, ϵ::SymmetricTensor, old::GeneralizedMaxwellState, Δt)
    σ0 = calculate_stress(m.base, ϵ)
    return mapreduce((c, ϵv) -> calculate_stress(c, ϵ, ϵv, Δt), +, m.chains, old.ϵv; init = σ0)
end

function MMB.material_response(m::GeneralizedMaxwell, ϵ::SymmetricTensor{2,3}, old::GeneralizedMaxwellState, Δt, _, _)
    dσdϵ, σ = gradient(e -> calculate_stress(m, e, old, Δt), ϵ, :all)
    state = GeneralizedMaxwellState(map((c, ϵv_old) -> calculate_viscous_strain(c, ϵ, ϵv_old, Δt), m.chains, old.ϵv))
    return σ, dσdϵ, state
end
