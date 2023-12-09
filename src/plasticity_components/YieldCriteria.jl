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

effective_stress(::VonMises, σred) = vonmises(σred)

effective_stress_gradient(yc::VonMises, σred) = (3/2)*dev(σred)/effective_stress(yc, σred) # More efficient than using AD above

MMB.get_num_params(::VonMises) = 1
MMB.vector2material(v::AbstractVector, ::VonMises; offset=0) = VonMises(v[offset+1])
function MMB.material2vector!(v::Vector, yl::VonMises; offset=0)
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

function effective_stress(yc::DruckerPrager, σred::SymmetricTensor{2,3,T}) where T
    return sqrt(T(3)/2)*norm(dev(σred)) - yc.B*tr(σred)
end

MMB.get_num_params(::DruckerPrager) = 2
MMB.vector2material(v::AbstractVector, ::DruckerPrager; offset=0) = DruckerPrager(v[offset+1], v[offset+2])
function MMB.material2vector!(v::Vector, yl::DruckerPrager; offset=0)
    v[offset+1] = yl.Y0
    v[offset+2] = yl.B
    return v
end

Base.show(io::IO, ::MIME"text/plain", yl::DruckerPrager) = println(io, "DruckerPrager with Y0=", yl.Y0, ", and B=", yl.B)