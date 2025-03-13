"""
    RateIndependent()

The evolution of the plastic multiplier for a rate-dependent material is given by 
the so-called KKT loading/unloading conditions
```math
\\dot{\\lambda} \\geq 0, \\quad \\varPhi \\leq 0, \\quad \\dot{\\lambda}\\varPhi = 0
```
"""
struct RateIndependent end

MMB.get_num_params(::RateIndependent) = 0
MMB.tovector!(v, ::RateIndependent; kwargs...) = nothing
MMB.fromvector(v, ::RateIndependent; kwargs...) = RateIndependent()

yield_residual(::RateIndependent, Φ, args...) = Φ

function Base.show(io::IO, ::MIME"text/plain", ::RateIndependent)
    println(io, "Rate independent response")
end

abstract type OverstressFunction end

yield_residual(ratelaw::OverstressFunction, Φ, Δλ, Δt, σy) = Δλ - Δt*overstress_function(ratelaw, Φ, σy)

"""
    NortonOverstress(; tstar, nexp)

The norton overstress function is defined as 
```math
\\eta(\\varPhi, \\sigma_\\mathrm{y}) = \\frac{1}{t_*} \\left\\langle \\frac{\\varPhi}{\\sigma_\\mathrm{y}} \\right\\rangle^n
```
where the material parameters ``t_*`` (`tstar`) and ``n`` (`nexp`) represent the 
relaxation time and overstress sensitivty.  
"""
@kwdef struct NortonOverstress{TT,TN} <: OverstressFunction
    tstar::TT
    nexp::TN
end

overstress_function(ratelaw::NortonOverstress, Φ, σy) = (1/ratelaw.tstar) * macaulay(Φ/σy)^ratelaw.nexp

MMB.get_num_params(::NortonOverstress) = 2
MMB.fromvector(v::AbstractVector, ::NortonOverstress; offset=0) = NortonOverstress(v[offset+1], v[offset+2])

function MMB.tovector!(v::Vector, o::NortonOverstress; offset=0)
    v[offset+1] = o.tstar
    v[offset+2] = o.nexp
    return v
end

function Base.show(io::IO, ::MIME"text/plain", o::NortonOverstress)
    println(io, "Norton overstress with t*=$(o.tstar) and n=$(o.nexp)")
end