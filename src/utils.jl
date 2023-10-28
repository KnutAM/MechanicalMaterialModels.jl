maketuple_or_nothing(::Nothing) = nothing
maketuple_or_nothing(x::Tuple) = x
maketuple_or_nothing(x::Any) = (x,)

get_promoted_type(args...) = promote_type(map(typeof, args)...)

""" 
    function vonmises(σ::SymmetricTensor{2,3})

Calculate the von Mises effective stress for a 2nd order tensor
"""
function vonmises(σ::SymmetricTensor{2,3,T}) where T
    return sqrt(T(3)/2)*norm(dev(σ))
end

"""
    function macaulay(x)

Calculate the macaulay bracket of x, ``\\langle x \\rangle``
```math
\\langle x \\rangle = \\left\\lbrace \\begin{matrix} 0 & x\\leq 0 \\\\ x & x>0 \\end{matrix} \\right .
```
"""
macaulay(x::T) where {T} = x > zero(T) ? x : zero(T)

function vector_residual!(R::Function, r_vector::AbstractVector{T}, x_vector::AbstractVector{T}, x) where T
    # construct residuals with type T
    x_tensor = frommandel(Tensors.get_base(typeof(x)), x_vector)
    r_tensor = R(x_tensor)
    tomandel!(r_vector, r_tensor)
    return r_vector
end

shapeof(t::TT) where{TT<:AbstractTensor} = Tensors.get_base(TT)