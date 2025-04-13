
order_params(x::NamedTuple{K}) where K = K[1]<K[2] ? x : reverse(x)

"""
    convert_hooke_param(T::Symbol; p1, p2)

Convert the hooke (isotropic) parameters `p1` and `p2` to the parameter with symbol `T`,
e.g. `G = convert_hooke_param(:G; E = 210e3, ν = 0.3)`
"""
convert_hooke_param(T::Symbol; kwargs...) = convert_hooke_param(Val(T); kwargs...)

function convert_hooke_param(::Val{T}; kwargs...) where {T}
    nt = NamedTuple(kwargs)
    p = order_params(nt)
    return convert_hooke_param(Val(T), Val(keys(p)); p...)
end

# Implemented for E, ν, G, K
convert_hooke_param(::Val{:E}, ::Val{(:G, :K)}; G, K) = 9K * G / (3K + G)
convert_hooke_param(::Val{:E}, ::Val{(:K, :ν)}; K, ν) = 3K * (1 - 2ν)
convert_hooke_param(::Val{:E}, ::Val{(:G, :ν)}; G, ν) = 2G * (1 + ν)

convert_hooke_param(::Val{:ν}, ::Val{(:E, :G)}; E, G) = (E / 2G) - 1
convert_hooke_param(::Val{:ν}, ::Val{(:E, :K)}; E, K) = (3K - E) / 6K
convert_hooke_param(::Val{:ν}, ::Val{(:G, :K)}; G, K) = (3K - 2G) / (2 * (3K + G))

convert_hooke_param(::Val{:G}, ::Val{(:E, :K)}; E, K) = (3K * E) / (9K - E)
convert_hooke_param(::Val{:G}, ::Val{(:E, :ν)}; E, ν) = E / (2 * (1 + ν))
convert_hooke_param(::Val{:G}, ::Val{(:K, :ν)}; K, ν) = (3K * (1 - 2ν)) / (2 * (1 + ν))

convert_hooke_param(::Val{:K}, ::Val{(:E, :G)}; E, G) = E * G / (3 * (3G - E))
convert_hooke_param(::Val{:K}, ::Val{(:E, :ν)}; E, ν) = E / (3 * (1 - 2ν))
convert_hooke_param(::Val{:K}, ::Val{(:G, :ν)}; G, ν) = 2G * (1 + ν) / (3 * (1 - 2ν))
