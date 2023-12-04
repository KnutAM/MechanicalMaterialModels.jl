# Isotropic hardening
abstract type AbstractIsotropicHardening{T} end
# Default initial value for κ
get_initial_value(::AbstractIsotropicHardening{T}) where {T} = zero(T) 

""" 
    Voce(;Hiso, κ∞)

Exponentially saturating isotropic hardening

```math
\\kappa_i = g_{\\mathrm{iso},i}(\\lambda) = \\kappa_\\infty \\left[1 - \\mathrm{exp}\\left(\\frac{H_\\mathrm{iso}}{\\kappa_\\infty} \\lambda \\right)\\right]
```
or alternatively as differential equations
```math
\\dot{\\kappa_i} = \\dot{\\lambda} H_\\mathrm{iso} \\left[1 - \\frac{\\kappa_i}{\\kappa_\\infty}\\right]
```

# Arguments
- `Hiso`: Isotropic hardening modulus, ``H_\\mathrm{iso}``
- `κ∞`: Saturation hardening value, ``\\kappa_\\infty``
    
"""
@kwdef struct Voce{T} <: AbstractIsotropicHardening{T}
    Hiso::T     # Initial hardening modulus
    κ∞::T       # Saturation stress
end

function get_evolution(param::Voce, κ::Number)
    return param.Hiso * (one(κ) - κ / param.κ∞)
end

""" 
    Swift(; K, λ0, n)

Isotropic hardening by the Swift power law

```math
\\kappa_i = g_{\\mathrm{iso},i}(\\lambda) = K \\left[\\lambda_0 + \\lambda \\right]^n
```

# Arguments
- `K`: ``K``
- `λ0`: ``\\lambda_0``
- `n`: ``n``

"""
@kwdef struct Swift{T} <:AbstractIsotropicHardening{T}
    K::T
    λ0::T
    n::T 
end

function get_evolution(param::Swift, κ::Number)
    return param.K * param.n * (κ/param.K)^((param.n-1)/param.n)
end

get_initial_value(param::Swift) = param.K*(param.λ0^param.n)


## Conversion from vector to material
MMB.get_num_params(::Voce) = 2
MMB.get_num_params(::Swift) = 3

MMB.vector2material(v::AbstractVector, ::Voce; offset=0) = Voce(Hiso=v[offset+1], κ∞=v[offset+2])
MMB.vector2material(v::AbstractVector, ::Swift; offset=0) = Swift(K=v[offset+1], λ0=v[offset+2], n=v[offset+3])

function MMB.material2vector!(v::AbstractVector, ih::IH; offset=0) where{IH<:AbstractIsotropicHardening}
    for (i, key) in enumerate(fieldnames(IH))
        v[offset+i] = getfield(ih, key)
    end
    return v
end

# Show methods 
function Base.show(io::IO, ::MIME"text/plain", h::Voce)
    println(io, "Voce with H=$(h.Hiso) and κ∞=$(h.κ∞)")
end
function Base.show(io::IO, ::MIME"text/plain", h::Swift)
    println(io, "Swift with K=$(h.K), λ0=$(h.λ0), and n=$(h.n)")
end
