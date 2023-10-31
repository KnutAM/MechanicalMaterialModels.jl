# Temporary solution untill Julia#51906 is solved
# https://github.com/JuliaLang/julia/issues/51906

"""
    LinearElastic(C::SymmetricTensor{4,3})
    LinearElastic{:general}(C::SymmetricTensor{4,3})

Creates a general linear elastic material (`symmetry=none`) where 
``\\boldsymbol{\\sigma} = \\boldsymbol{C}:\\boldsymbol{\\epsilon}``. 
``\\boldsymbol{C}`` is the 4th order elastic stiffness tensor. 

# Arguments
- `C::SymmetricTensor{4,3}`: 4th order symmetric elastic stiffness tensor


    LinearElastic(; E, ν)
    LinearElastic{:isotropic}(; E, ν)

Creates an isotropic LinearElastic material, such that
```math
\\boldsymbol{\\sigma} = 2\\mu \\boldsymbol{\\epsilon} + \\lambda \\mathrm{tr}(\\boldsymbol{\\epsilon}) \\boldsymbol{I} \\\\
```
where 
```math
\\mu = \\frac{E}{2(1+\\nu)}, \\quad \\lambda=\\frac{E\\nu}{(1+\\nu)(1-2\\nu)}
```

# Arguments
- `E`: Young's modulus, ``E``
- `ν`: Poisson's ratio, ``\\nu``


    LinearElastic{:cubicsymmetry}(; C1111::T, C1122::T, C1212::T) where {T}

Creates a LinearElastic type, hence ``\\boldsymbol{\\sigma} = \\boldsymbol{C}:\\boldsymbol{\\epsilon}``, where
``\\boldsymbol{C}`` possesses cubic symmetry along the coordinate axes. It is defined as

```math
\\boldsymbol{C} = 
\\begin{bmatrix}
C_{1111} & C_{1122} & C_{1122} & 0 & 0 & 0 & 0 & 0 & 0 \\\\
C_{1122} & C_{1111} & C_{1122} & 0 & 0 & 0 & 0 & 0 & 0 \\\\
C_{1122} & C_{1122} & C_{1111} & 0 & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 2C_{1212} & 0 & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 2C_{1212} & 0 & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 2C_{1212} & 0 & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 2C_{1212} & 0 & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 2C_{1212} & 0 \\\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2C_{1212} \\\\
\\end{bmatrix}
```

# Arguments
- `C1111`: The ``C_{1111}`` stiffness parameter
- `C1122`: The ``C_{1122}`` stiffness parameter
- `C1212`: The ``C_{1212}`` stiffness parameter
"""
LinearElastic

struct LinearElastic{T, case, N} <: AbstractMaterial
    C::SymmetricTensor{4,3,T,36}
    p::SVector{N,T}
end

get_base_numbertype(::LinearElastic{T}) where T = T

# General symmetry
LinearElastic{:general}(C::SymmetricTensor) = LinearElastic(C)
LinearElastic(C::SymmetricTensor{4,3,T}) where{T} = LinearElastic{T,:general,36}(C,SVector{36,T}(reshape(tomandel(C),36)))

# Isotropic
function LinearElastic{:isotropic}(;E, ν)
    T = get_promoted_type(E, ν)
    μ = E/(2*(1+ν))
    λ = E*ν/((1+ν)*(1-2*ν))
    I2 = one(SymmetricTensor{2,3,T})
    C = symmetric(2μ * otimesu(I2,I2) + λ * I2⊗I2)
    return LinearElastic{T, :isotropic, 2}(C, SVector(T(E), T(ν)))
end
LinearElastic(; E, ν) = LinearElastic{:isotropic}(; E, ν)

# Cubic symmetry
function LinearElastic{:cubicsymmetry}(; C1111, C1122, C1212)
    T = get_promoted_type(C1111, C1122, C1212)
    C = SymmetricTensor{4,3}(
        function f(i,j,k,l)
            i==j==k==l     && return C1111
            (i==j && k==l) && return C1122
            (i,j)==(k,l)   && return 2*C1212
            return zero(T)
        end
        )
    return LinearElastic{T, :cubicsymmetry, 3}(C, SVector(C1111, C1122, C1212))
end

function MMB.material_response(m::LinearElastic, ϵ::SymmetricTensor{2,3}, old=NoMaterialState(), Δt=nothing, cache=allocate_material_cache(m), extras::AbstractExtraOutput=NoExtraOutput(); options=Dict{Symbol, Any}())
    σ = calculate_stress(m, ϵ)
    return σ, m.C, old
end

calculate_stress(m::LinearElastic, ϵ::SymmetricTensor) = m.C⊡ϵ

# Functions for conversion between material and parameter vectors
MMB.get_num_params(::LinearElastic{<:Any,<:Any,N}) where{N} = N

function MMB.vector2material(v::AbstractVector, ::LinearElastic{<:Any, :isotropic}; offset=0)
    return LinearElastic{:isotropic}(;E=v[offset+1], ν=v[offset+2])
end
function MMB.vector2material(v::AbstractVector, ::LinearElastic{<:Any, :cubicsymmetry}; offset=0)
    return LinearElastic{:cubicsymmetry}(;C1111=v[offset+1], C1122=v[offset+2], C1212=v[offset+3])
end

function MMB.material2vector!(v::AbstractVector, m::LinearElastic; offset=0)
    for (i, p) in pairs(m.p)
        v[offset+i] = p
    end
    return v
end

# Differentiation
function MMB.differentiate_material!(deriv::MaterialDerivatives{T}, m::LinearElastic, ϵ, ⁿs, Δt, dσdϵ, ::AbstractExtraOutput, ::Any) where {T}
    tomandel!(deriv.dσdϵ, dσdϵ)
    p = material2vector(m)
    σ_from_param(p_vector) = tomandel(calculate_stress(vector2material(p_vector, m), ϵ))
    ForwardDiff.jacobian!(deriv.dσdp,σ_from_param,p)
    
    fill!(deriv.dσdⁿs,  zero(T))
    fill!(deriv.dsdϵ,  zero(T))
    fill!(deriv.dsdp,  zero(T))
    fill!(deriv.dsdⁿs,  zero(T))
end

# Show methods
function Base.show(io::IO, M::MIME"text/plain", e::LinearElastic{<:Any,case}) where case
    print(io, "LinearElastic: ")
    if case === :isotropic
        print(io, "Isotropic with E=$(e.p[1]) and ν=$(e.p[2])")
    elseif case === :cubicsymmetry
        print(io, "Cubic symmetric with  C1111=$(e.p[1]), C1122=$(e.p[2]), C1212=$(e.p[3])")
    elseif case === :general
        print(io, "Fully anisotropic with mandel stiffness tensor")
        show(io, M, tomandel(e.C))
    else
        print(io, "$case with mandel stiffness tensor")
        show(io, M, tomandel(e.C))
    end
    println(io)
end
