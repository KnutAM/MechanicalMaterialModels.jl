"""
    maketuple_or_nothing(x)

* `x` is a single value: Convert to a `Tuple` of length 1
* `x` is a `Tuple` or `Nothing`: Return `x`
"""
maketuple_or_nothing(::Nothing) = nothing
maketuple_or_nothing(x::Tuple) = x
maketuple_or_nothing(x::Any) = (x,)

"""
    get_promoted_type(args...)

Get the promoted type for the type of the arguments, 
e.g. `get_promoted_type(1, 1.f0)` is `Float32`
"""
get_promoted_type(args...) = promote_type(map(typeof, args)...)

""" 
    function vonmises(σ::SecondOrderTensor{3})

Calculate the von Mises effective stress for a 2nd order tensor
"""
function vonmises(σ::SecondOrderTensor{3})
    σ_dev = dev(σ)
    return sqrt(3*(σ_dev ⊡ transpose(σ_dev))/2)
end

""" 
    function vonmises_and_gradient(σ::SecondOrderTensor{3})

Calculate the von Mises effective stress for a 2nd order tensor
as well as the gradient, ν
"""
function vonmises_and_gradient(σ::SecondOrderTensor{3})
    σ_dev = dev(σ)
    σvm = sqrt(3*(σ_dev ⊡ transpose(σ_dev))/2)
    return σvm, (3/(2*σvm))*transpose(σ_dev)
end

"""
    function macaulay(x)

Calculate the macaulay bracket of x, ``\\langle x \\rangle``
```math
\\langle x \\rangle = \\left\\lbrace \\begin{matrix} 0 & x\\leq 0 \\\\ x & x>0 \\end{matrix} \\right .
```
"""
macaulay(x::T) where {T} = x > zero(T) ? x : zero(T)

"""
    baseof(t::AbstractTensor)

Get the base-type of `t`, i.e. if `t::SymmetricTensor{2,3,Float64,6}`,
`baseof(t)` returns `SymmetricTensor{2,3}`
"""
baseof(::TT) where{TT} = Tensors.get_base(TT)

"""
    vector_residual!(rf::Function, r_vector::AbstractVector, x_vector::AbstractVector, x)

Makes it easy to construct a mutating vector residual function from a tensor-like equation,
`r = rf(x) = residual(x, args...)`, e.g.
`rf!(r_vector, x_vector) = vector_residual!(z -> residual(z, args...), r_vector, x_vector, x)`

The input `x` and output `r` of `rf` should have the same type, `RT`, and support
`Tensors.get_base(RT)`, `frommandel(Tensors.get_base(RT), x_vector)`, and 
`tomandel!(r_vector, r)`.

The approach was adopted from https://github.com/kimauth/MaterialModels.jl
"""
function vector_residual!(rf::F, r_vector::AbstractVector{T}, x_vector::AbstractVector{T}, x) where {F<:Function, T}
    x_tensor = frommandel(baseof(x), x_vector)
    r_tensor = rf(x_tensor)
    tomandel!(r_vector, r_tensor)
    return r_vector
end

"""
    static_vector(args::Union{Number, SVector}...)

Convert the elements of `args` into an `SVector`.
"""
@inline static_vector(a, b, c...) = static_vector(static_vector(a,b), static_vector(c...))
@inline static_vector(a::Number, b::Number) = static_vector(SVector(a), SVector(b))
@inline static_vector(a::Number, b::StaticVector) = static_vector(SVector(a), b)
@inline static_vector(a::StaticVector, b::Number) = static_vector(a, SVector(b))

@inline static_vector(a::StaticVector, b::StaticVector) = vcat(a, b)
@inline static_vector(a::StaticVector) = a