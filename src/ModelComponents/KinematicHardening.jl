# Kinematic hardening
abstract type AbstractKinematicHardening{T} end

"""
    ArmstrongFrederick(Hkin, β∞)

Armstrong-Frederick kinematic hardening law (doi: 10.1179/096034007X207589)

```math
g_{\\mathrm{kin},i}(\\nu, \\boldsymbol{\\beta}_i) = H_\\mathrm{kin} (\\frac{2}{3}\\boldsymbol{\\nu} - \\frac{\\boldsymbol{\\beta}_i}{\\beta_\\infty})
```

# Arguments
- `Hkin`: Kinematic hardening modulus, ``H_\\mathrm{kin}``
- `β∞`: Effective back-stress saturation value, ``\\beta_\\infty``
"""
struct ArmstrongFrederick{T} <: AbstractKinematicHardening{T}
    Hkin::T     # Initial hardening modulus
    β∞::T       # Saturation stress
end
ArmstrongFrederick(;Hkin, β∞) = ArmstrongFrederick(Hkin, β∞)    # Keyword argument constructor

function get_evolution(param::ArmstrongFrederick, ν::SecondOrderTensor, βᵢ::SecondOrderTensor)
    # (2/3)*Hkin(ν - (3/2)*βᵢ/β∞)
    param.Hkin * ((2.0/3.0) * ν - βᵢ/param.β∞)
end

"""
    Delobelle(Hkin, β∞, δ)

Kinematic hardening law according to Delobelle, which combines the Armstrong-Frederick law with the Burlet-Cailletaud law
(doi: 10.1016/S0749-6419(95)00001-1)

```math
g_{\\mathrm{kin},i}(\\nu, \\boldsymbol{\\beta}_i) = H_\\mathrm{kin} \\left[\\frac{2}{3}\\boldsymbol{\\nu} 
                                    - \\delta\\frac{\\boldsymbol{\\beta}_i}{\\beta_\\infty}
                                    - \\frac{2}{3\\beta_\\infty}\\left[1 - \\delta\\right]\\left[\\boldsymbol{\\nu}:\\boldsymbol{\\beta}_i\\right]\\boldsymbol{\\nu}
                                    \\right]
```

# Arguments
- `Hkin`: Kinematic hardening modulus, ``H_\\mathrm{kin}``
- `β∞`: Effective back-stress saturation value, ``\\beta_\\infty``
- `δ`: Amount of Armstrong-Frederick type of kinematic hardening, ``\\delta``

"""
struct Delobelle{T} <: AbstractKinematicHardening{T}
    Hkin::T     # Initial hardening modulus
    β∞::T       # Saturation stress
    δ::T        # Amount of Armstrong-Frederick hardening
end
Delobelle(;Hkin, β∞, δ) = Delobelle(Hkin, β∞, δ)    # Keyword argument constructor

function get_evolution(param::Delobelle, ν::SecondOrderTensor, βᵢ::SecondOrderTensor)
    AF_Term = (param.δ/param.β∞) * βᵢ                        # Armstrong Frederick term
    BC_Term = (2.0/3.0) * (1.0-param.δ)*((ν⊡βᵢ)/param.β∞)*ν  # Burlet Cailletaud term
    return param.Hkin * ((2.0/3.0) * ν - AF_Term - BC_Term)  # Complete evolution 
end


"""
    OhnoWang(Hkin, β∞, m)

Kinematic hardening law according to Ohno-Wang (doi: 10.1016/0749-6419(93)90042-O)

```math
g_{\\mathrm{kin},i}(\\nu, \\boldsymbol{\\beta}_i) = H_\\mathrm{kin} \\left[\\frac{2}{3}\\boldsymbol{\\nu} 
                                    - \\frac{\\boldsymbol{\\beta}_i}{\\beta_\\infty} 
                                    \\frac{\\langle \\boldsymbol{\\nu}:\\boldsymbol{\\beta}_i \\rangle}{\\beta_\\infty}
                                    \\left[\\frac{\\beta_i^\\mathrm{vM}}{\\beta_\\infty}\\right]^m
                                    \\right]
```
where ``\\langle x \\rangle`` is 0 if ``x\\leq 0`` and ``x`` if ``x>0``.
``\\beta_i^\\mathrm{vM} = \\sqrt{2\\boldsymbol{\\beta}_i:\\boldsymbol{\\beta}_i/3}``, noting that
``\\boldsymbol{\\beta}_i`` is deviatoric.

# Arguments
- `Hkin`: Kinematic hardening modulus, ``H_\\mathrm{kin}``
- `β∞`: Effective back-stress saturation value, ``\\beta_\\infty``
- `m`: Exponent in the OhnoWang equation, ``m``

"""
struct OhnoWang{T} <: AbstractKinematicHardening{T}
    Hkin::T     # Initial hardening modulus
    β∞::T       # Saturation stress
    m::T     # Ohno Wang exponent
end
OhnoWang(;Hkin, β∞, m) = OhnoWang(Hkin, β∞, m)    # Keyword argument constructor

function get_evolution(param::OhnoWang{Tp}, ν::SecondOrderTensor, βᵢ::SecondOrderTensor{dim,Tβ}) where{Tp,Tβ,dim}
    β_vm = vonmises(βᵢ)
    if β_vm < param.β∞ * eps(promote_type(Tp,Tβ))
        return param.Hkin * (2.0/3.0) * ν + 0*βᵢ
    end
    mac_term = (macaulay(ν⊡βᵢ) /param.β∞)
    exp_term = (β_vm/param.β∞)^param.m
    return param.Hkin * ((2.0/3.0) * ν - βᵢ * mac_term * exp_term / β_vm )
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