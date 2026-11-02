abstract type YieldCriterion end

function yield_criterion(yc::YieldCriterion, σred::AbstractTensor, ΔY::Number)
    return effective_stress(yc, σred) - (initial_yield_limit(yc)+ΔY)
end

effective_stress_gradient(yc, σred) = gradient(σ -> effective_stress(yc, σ), σred)

default_yield_criteria(yc::YieldCriterion) = yc
default_yield_criteria(v::Real) = VonMises(v)

"""
    VonMises(; Y0)

Create a von Mises yield criterion with initial yield limit, ``Y_0``, as `Y0`.
The yield criterion is then defined as
```math
\\Phi = \\sqrt{\\frac{3}{2}} \\left| \\text{dev} \\left( \\boldsymbol{\\sigma}_\\mathrm{red} \\right) \\right| - \\left[ Y_0 + \\Delta Y \\right] = 0
```
where ``\\boldsymbol{\\sigma}_\\mathrm{red}`` is the reduced (by kinematic hardening) stress tensor, and ``\\Delta Y`` the change of the initial 
yield limit due to isotropic hardening (i.e. ``\\kappa``).
"""
@kwdef struct VonMises{T} <: YieldCriterion
    Y0::T
end
initial_yield_limit(yc::VonMises) = yc.Y0

effective_stress(::VonMises, σred::SecondOrderTensor{3}) = vonmises(σred)

 # More efficient than using AD above
function effective_stress_gradient(::VonMises, σred::SecondOrderTensor{3})
    σred_dev = dev(σred)
    σred_dev_t = transpose(σred_dev)
    return σred_dev_t / sqrt(2*(σred_dev ⊡ σred_dev_t)/3)
end

MMB.get_vector_length(::VonMises) = 1
MMB.get_vector_eltype(::VonMises{T}) where {T} = T
MMB.fromvector(v::AbstractVector, ::VonMises; offset=0) = VonMises(v[offset+1])
function MMB.tovector!(v::Vector, yl::VonMises; offset=0)
    v[offset+1] = yl.Y0
    return v
end

Base.show(io::IO, ::MIME"text/plain", yl::VonMises) = println(io, "VonMises with Y0=", yl.Y0)

"""
    DruckerPrager(; Y0, B)

Create a Drucker-Prager yield criterion, with initial yield limit, ``Y_0``, as `Y0`,
and pressure sensitivity `B`. The yield criterion is defined as 
```math
\\Phi = \\sqrt{\\frac{3}{2}} \\left| \\mathrm{dev} \\left( \\boldsymbol{\\sigma}_\\mathrm{red} \\right) \\right| - B\\mathrm{tr}\\left( \\boldsymbol{\\sigma}_\\mathrm{red} \\right) - \\left[ Y_0 + \\Delta Y \\right] = 0
```
where ``\\boldsymbol{\\sigma}_\\mathrm{red}`` is the reduced (by kinematic hardening) stress tensor, and ``\\Delta Y`` the change of the initial 
yield limit due to isotropic hardening (i.e. ``\\kappa``).
"""
@kwdef struct DruckerPrager{T} <: YieldCriterion
    Y0::T
    B::T
end

initial_yield_limit(yc::DruckerPrager) = yc.Y0

function effective_stress(yc::DruckerPrager, σred::SecondOrderTensor{3})
    return vonmises(σred) - yc.B*tr(σred)
end

# More efficient than using AD above, and handles sub-gradient at dev(σred) = 0
function effective_stress_gradient(yc::DruckerPrager, σred::SecondOrderTensor{3})
    σred_dev = dev(σred)
    σred_dev_t = transpose(σred_dev)
    σred_dev_sq = σred_dev ⊡ σred_dev_t
    ν_tr = - yc.B * one(σred)
    iszero(σred_dev_sq) && return ν_tr # Subdifferential
    ν_dev = σred_dev_t / sqrt(2 * σred_dev_sq / 3)
    return ν_dev + ν_tr
end

MMB.get_vector_length(::DruckerPrager) = 2
MMB.get_vector_eltype(::DruckerPrager{T}) where {T} = T
MMB.fromvector(v::AbstractVector, ::DruckerPrager; offset=0) = DruckerPrager(v[offset+1], v[offset+2])
function MMB.tovector!(v::Vector, yl::DruckerPrager; offset=0)
    v[offset+1] = yl.Y0
    v[offset+2] = yl.B
    return v
end

Base.show(io::IO, ::MIME"text/plain", yl::DruckerPrager) = println(io, "DruckerPrager with Y0=", yl.Y0, ", and B=", yl.B)