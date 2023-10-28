# The evolution of the plastic multiplier for a rate-dependent material is 
# ∂λ/∂t = f(Φ)
# With backward Euler, 
# λ-λⁿ = Δt*f(Φ) with residual 
# rΦ = λ- (λⁿ+Δt*f(Φ))
struct RateIndependent end

MMB.get_num_params(::RateIndependent) = 0
MMB.material2vector!(v, ::RateIndependent; kwargs...) = nothing
MMB.vector2material(v, ::RateIndependent; kwargs...) = RateIndependent()

yield_residual(::RateIndependent, Φ, args...) = Φ

function Base.show(io::IO, ::MIME"text/plain", ::RateIndependent)
    println(io, "Rate independent response")
end

abstract type OverstressFunction end

yield_residual(ratelaw::OverstressFunction, Φ, Δλ, Δt, σy) = Δλ - Δt*overstress_function(ratelaw, Φ, σy)

struct NortonOverstress{TT,TN} <: OverstressFunction
    tstar::TT
    nexp::TN
end
NortonOverstress(;tstar, nexp) = NortonOverstress(tstar, nexp)

overstress_function(ratelaw::NortonOverstress, Φ, σy) = Φ<=0 ? zero(Φ) : (1/ratelaw.tstar)*(Φ/σy)^ratelaw.nexp

MMB.get_num_params(::NortonOverstress) = 2
MMB.vector2material(v::AbstractVector, ::NortonOverstress; offset=0) = NortonOverstress(v[offset+1], v[offset+2])

function MMB.material2vector!(v::Vector, o::NortonOverstress; offset=0)
    v[offset+1] = o.tstar
    v[offset+2] = o.nexp
    return v
end

function Base.show(io::IO, ::MIME"text/plain", o::NortonOverstress)
    println(io, "Norton overstress with t*=$(o.tstar) and n=$(o.nexp)")
end