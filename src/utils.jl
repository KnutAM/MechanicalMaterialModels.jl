"""
    get_promoted_type(args...)

Get the promoted type for the type of the arguments,
e.g. `get_promoted_type(1, 1.f0)` is `Float32`
"""
get_promoted_type(args...) = promote_type(map(typeof, args)...)

"""
    promote_tuple_eltype(::Type{T}, xs::Tuple) where {T}

Fold over the entries of `xs` (`Number`s or `Tensors.jl` tensors), promoting `T`
with each entry's scalar element type. Returns `T` unchanged for an empty tuple.
"""
promote_tuple_eltype(::Type{T}, ::Tuple{}) where {T} = T
function promote_tuple_eltype(::Type{T}, xs::Tuple) where {T}
    return promote_tuple_eltype(promote_type(T, eltype(first(xs))), Base.tail(xs))
end

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
    σ_dev_t = transpose(σ_dev)
    σvm = sqrt(3*(σ_dev ⊡ σ_dev_t)/2)
    return σvm, (3 / (2 * σvm)) * σ_dev_t
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

abstract type AbstractResidual end

"""
    vector_residual!(rf::Function, r_vector::AbstractVector, x_vector::AbstractVector, xbase::RT)

Makes it easy to construct a mutating vector residual function from a tensor-like equation,
`r = rf(x) = residual(x, args...)`, e.g.
`rf!(r_vector, x_vector) = vector_residual!(z -> residual(z, args...), r_vector, x_vector, x)`

The input `x::RT` must support MaterialModelsBase.fromvector, and the output `r` from `rf` 
must support`MaterialModelsBase.fromvector!`.

The approach was adopted from [MaterialModels.jl](https://github.com/kimauth/MaterialModels.jl)
"""
function vector_residual!(rf::F, r_vector::AbstractVector{T}, x_vector::AbstractVector{T}, xbase) where {F<:Function, T}
    x_tensor = fromvector(x_vector, xbase)
    r_tensor = rf(x_tensor)
    tovector!(r_vector, r_tensor)
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

"""
    remove_dual(x::T) where {T <: Number}
    remove_dual(x::AbstractTensor{<:Any, <:Any, T}) where {T}

Removes the dual part if `T <: ForwardDiff.Dual`, extract the value part, helpful
when implementing custom gradients.
"""
function remove_dual end

# Scalars
remove_dual(x::ForwardDiff.Dual) = ForwardDiff.value(x)
remove_dual(x::Number) = x
# Tensors
remove_dual(x::AbstractTensor{<:Any, <:Any, <:ForwardDiff.Dual}) = Tensors._extract_value(x)
remove_dual(x::AbstractTensor) = x
