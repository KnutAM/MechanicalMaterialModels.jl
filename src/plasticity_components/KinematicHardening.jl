# Kinematic hardening
abstract type AbstractKinematicHardening{T} end

"""
    ArmstrongFrederick(; Hkin, β∞)

Armstrong-Frederick kinematic hardening law with modulus `Hkin` and saturation
stress `β∞`. 
"""
@kwdef struct ArmstrongFrederick{T} <: AbstractKinematicHardening{T}
    Hkin::T     # Initial hardening modulus
    β∞::T       # Saturation stress
end

get_modulus(h::ArmstrongFrederick) = h.Hkin

function get_evolution(param::ArmstrongFrederick, ν::SecondOrderTensor, βᵢ::SecondOrderTensor)
    return - ν + transpose(βᵢ)*(3/(2*param.β∞))
end

"""
    Delobelle(; Hkin, β∞, δ)

Kinematic hardening law according to Delobelle with hardening modulus `Hkin`,
saturation stress, `β∞`, and scaling parameter `δ`, which scales between pure
Armstrong-Frederick hardening, `δ=1`, and Burlet-Cailletaud hardening, `δ=0`.
"""
@kwdef struct Delobelle{T} <: AbstractKinematicHardening{T}
    Hkin::T     # Initial hardening modulus
    β∞::T       # Saturation stress
    δ::T        # Amount of Armstrong-Frederick hardening
end

get_modulus(h::Delobelle) = h.Hkin

function get_evolution(param::Delobelle, ν::SecondOrderTensor, βᵢ::SecondOrderTensor)
    AF_Term = (3 * param.δ / (2 * param.β∞)) * transpose(βᵢ) # Armstrong Frederick saturation
    BC_Term = (1 - param.δ) * ((ν ⊡ βᵢ) / param.β∞) * ν     # Burlet Cailletaud saturation
    return - ν + AF_Term + BC_Term
end


"""
    OhnoWang(; Hkin, β∞, m)

Kinematic hardening law according to Ohno-Wang with hardening `Hkin`,
saturation stress, `β∞`, and exponent, `m`.
"""
@kwdef struct OhnoWang{T} <: AbstractKinematicHardening{T}
    Hkin::T  # Initial hardening modulus
    β∞::T    # Saturation stress
    m::T     # Ohno Wang exponent
end

get_modulus(h::OhnoWang) = h.Hkin

function get_evolution(param::OhnoWang{Tp}, ν::SecondOrderTensor, βᵢ::SecondOrderTensor{dim,Tβ}) where{Tp,Tβ,dim}
    β_vm = vonmises(βᵢ)
    if β_vm < param.β∞ * eps(promote_type(Tp,Tβ))
        return - ν + 0*βᵢ
    end
    mac_term = macaulay(3 * (ν ⊡ βᵢ)) / (2*param.β∞)
    exp_term = (β_vm / param.β∞)^param.m
    return - ν + transpose(βᵢ) * (mac_term * exp_term / β_vm )
end

## Conversion from vector to material
MMB.get_num_params(::ArmstrongFrederick) = 2
MMB.get_num_params(::Union{Delobelle,OhnoWang}) = 3

MMB.vector2material(v::AbstractVector, ::ArmstrongFrederick; offset=0) = ArmstrongFrederick(Hkin=v[offset+1], β∞=v[offset+2])
MMB.vector2material(v::AbstractVector, ::Delobelle; offset=0) = Delobelle(Hkin=v[offset+1], β∞=v[offset+2], δ=v[offset+3])
MMB.vector2material(v::AbstractVector, ::OhnoWang; offset=0) = OhnoWang(Hkin=v[offset+1], β∞=v[offset+2], m=v[offset+3])

function MMB.material2vector!(v::Vector, kh::KH; offset=0) where{KH<:AbstractKinematicHardening}
    for (i,key) in enumerate(fieldnames(KH))
        v[offset+i] = getfield(kh, key)
    end
    return v
end

# Show methods 
function Base.show(io::IO, M::MIME"text/plain", h::ArmstrongFrederick)
    println(io, "ArmstrongFrederick with H=$(h.Hkin) and β∞=$(h.β∞)")
end
function Base.show(io::IO, M::MIME"text/plain", h::Delobelle)
    println(io, "Delobelle with H=$(h.Hkin), β∞=$(h.β∞), and δ=$(h.δ)")
end
function Base.show(io::IO, M::MIME"text/plain", h::OhnoWang)
    println(io, "OhnoWang with H=$(h.Hkin), β∞=$(h.β∞), and m=$(h.m)")
end