"Get the slip systems for the given crystallography type."
function get_slip_systems end

# Construction of a `Crystallography` subtype, `CT`, are assumed support `CT(planes::SVector, directions::SVector)`
abstract type Crystallography{T} end
(::Type{CT})() where {CT<:Crystallography} = CT{Float64}()  # Default to Float64 slip definitions
(::Type{CT})() where {CT<:Crystallography{T}} where {T<:AbstractFloat} = CT(get_slip_systems(CT)...)

get_num_slipsystems(c::Crystallography) = length(get_slip_planes(c))
get_slip_plane(c::Crystallography, i::Int) = get_slip_planes(c)[i]
get_slip_direction(c::Crystallography, i::Int) = get_slip_directions(c)[i]
get_slip_dyads(c::Crystallography) = map(⊗, get_slip_directions(c), get_slip_planes(c))
get_slip_dyads_test(c::Crystallography) = map((s,m) -> (s⊗m), get_slip_directions(c), get_slip_planes(c))

struct BCC{T} <: Crystallography{T}
    planes::SVector{48,Vec{3,T}}
    directions::SVector{48,Vec{3,T}}
end

struct FCC{T} <: Crystallography{T}
    planes::SVector{12,Vec{3,T}}
    directions::SVector{12,Vec{3,T}}
end

struct BCC12{T} <: Crystallography{T}
    planes::SVector{12,Vec{3,T}}
    directions::SVector{12,Vec{3,T}}
end

struct GenericCrystallography{T,dim,N} <: Crystallography{T}
    planes::SVector{N, Vec{dim,T}}
    directions::SVector{N, Vec{dim,T}}
end

GenericCrystallography(args...) = GenericCrystallography{Float64}(args...)

function GenericCrystallography{T}(angles::Number...) where {T<:AbstractFloat}
    p = Vec{2,T}((0, 1))
    d = Vec{2,T}((1, 0))
    return GenericCrystallography(
        SVector(map(θ -> rotate(p, θ), angles)),
        SVector(map(θ -> rotate(d, θ), angles)))
end
# TODO: Extend to 3d when `eltype(angles)<:Tensor{2,3}`

get_slip_planes(c::Union{BCC,FCC,BCC12,GenericCrystallography}) = c.planes
get_slip_directions(c::Union{BCC,FCC,BCC12,GenericCrystallography}) = c.directions

# Generated with the file "python/slip_systems.py",
# adapted from from github.com/knutam/MaterialModels
function get_slip_systems(::Type{BCC{T}}) where T
    slip_planes = SVector{48,Vec{3,T}}(
        Vec{3,T}(( 1,  1,  0))/sqrt(2),    # 1
        Vec{3,T}(( 1,  1,  0))/sqrt(2),    # 2
        Vec{3,T}(( 1, -1,  0))/sqrt(2),    # 3
        Vec{3,T}(( 1, -1,  0))/sqrt(2),    # 4
        Vec{3,T}(( 1,  0,  1))/sqrt(2),    # 5
        Vec{3,T}(( 1,  0,  1))/sqrt(2),    # 6
        Vec{3,T}(( 1,  0, -1))/sqrt(2),    # 7
        Vec{3,T}(( 1,  0, -1))/sqrt(2),    # 8
        Vec{3,T}(( 0,  1,  1))/sqrt(2),    # 9
        Vec{3,T}(( 0,  1,  1))/sqrt(2),    # 10
        Vec{3,T}(( 0,  1, -1))/sqrt(2),    # 11
        Vec{3,T}(( 0,  1, -1))/sqrt(2),    # 12
        Vec{3,T}(( 1,  1,  2))/sqrt(6),    # 13
        Vec{3,T}(( 1,  1, -2))/sqrt(6),    # 14
        Vec{3,T}(( 1, -1,  2))/sqrt(6),    # 15
        Vec{3,T}(( 1, -1, -2))/sqrt(6),    # 16
        Vec{3,T}(( 1,  2,  1))/sqrt(6),    # 17
        Vec{3,T}(( 1,  2, -1))/sqrt(6),    # 18
        Vec{3,T}(( 1, -2,  1))/sqrt(6),    # 19
        Vec{3,T}(( 1, -2, -1))/sqrt(6),    # 20
        Vec{3,T}(( 2,  1,  1))/sqrt(6),    # 21
        Vec{3,T}(( 2,  1, -1))/sqrt(6),    # 22
        Vec{3,T}(( 2, -1,  1))/sqrt(6),    # 23
        Vec{3,T}(( 2, -1, -1))/sqrt(6),    # 24
        Vec{3,T}(( 1,  2,  3))/sqrt(14),   # 25
        Vec{3,T}(( 1,  2, -3))/sqrt(14),   # 26
        Vec{3,T}(( 1, -2,  3))/sqrt(14),   # 27
        Vec{3,T}(( 1, -2, -3))/sqrt(14),   # 28
        Vec{3,T}(( 1,  3,  2))/sqrt(14),   # 29
        Vec{3,T}(( 1,  3, -2))/sqrt(14),   # 30
        Vec{3,T}(( 1, -3,  2))/sqrt(14),   # 31
        Vec{3,T}(( 1, -3, -2))/sqrt(14),   # 32
        Vec{3,T}(( 2,  1,  3))/sqrt(14),   # 33
        Vec{3,T}(( 2,  1, -3))/sqrt(14),   # 34
        Vec{3,T}(( 2, -1,  3))/sqrt(14),   # 35
        Vec{3,T}(( 2, -1, -3))/sqrt(14),   # 36
        Vec{3,T}(( 2,  3,  1))/sqrt(14),   # 37
        Vec{3,T}(( 2,  3, -1))/sqrt(14),   # 38
        Vec{3,T}(( 2, -3,  1))/sqrt(14),   # 39
        Vec{3,T}(( 2, -3, -1))/sqrt(14),   # 40
        Vec{3,T}(( 3,  1,  2))/sqrt(14),   # 41
        Vec{3,T}(( 3,  1, -2))/sqrt(14),   # 42
        Vec{3,T}(( 3, -1,  2))/sqrt(14),   # 43
        Vec{3,T}(( 3, -1, -2))/sqrt(14),   # 44
        Vec{3,T}(( 3,  2,  1))/sqrt(14),   # 45
        Vec{3,T}(( 3,  2, -1))/sqrt(14),   # 46
        Vec{3,T}(( 3, -2,  1))/sqrt(14),   # 47
        Vec{3,T}(( 3, -2, -1))/sqrt(14),   # 48
        )
    slip_directions = SVector{48,Vec{3,T}}(
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 1
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 2
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 3
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 4
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 5
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 6
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 7
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 8
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 9
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 10
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 11
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 12
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 13
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 14
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 15
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 16
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 17
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 18
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 19
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 20
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 21
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 22
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 23
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 24
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 25
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 26
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 27
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 28
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 29
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 30
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 31
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 32
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 33
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 34
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 35
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 36
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 37
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 38
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 39
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 40
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 41
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 42
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 43
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 44
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 45
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 46
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 47
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 48
    )
    return slip_planes, slip_directions
end

function get_slip_systems(::Type{FCC{T}}) where T
    slip_planes = SVector{12,Vec{3,T}}(
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 1 
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 2 
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 3 
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 4 
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 5 
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 6 
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 7 
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 8 
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 9 
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 10
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 11
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 12
        )
    slip_directions = SVector{12,Vec{3,T}}(
        Vec{3,T}(( 1, -1,  0))/sqrt(2),    # 1 
        Vec{3,T}(( 1,  0, -1))/sqrt(2),    # 2 
        Vec{3,T}(( 0,  1, -1))/sqrt(2),    # 3 
        Vec{3,T}(( 1, -1,  0))/sqrt(2),    # 4 
        Vec{3,T}(( 1,  0,  1))/sqrt(2),    # 5 
        Vec{3,T}(( 0,  1,  1))/sqrt(2),    # 6 
        Vec{3,T}(( 1,  1,  0))/sqrt(2),    # 7 
        Vec{3,T}(( 1,  0, -1))/sqrt(2),    # 8 
        Vec{3,T}(( 0,  1,  1))/sqrt(2),    # 9 
        Vec{3,T}(( 1,  1,  0))/sqrt(2),    # 10
        Vec{3,T}(( 1,  0,  1))/sqrt(2),    # 11
        Vec{3,T}(( 0,  1, -1))/sqrt(2),    # 12
    )
    return slip_planes, slip_directions
end

function get_slip_systems(::Type{BCC12{T}}) where T
    slip_planes = SVector{12,Vec{3,T}}(
        Vec{3,T}(( 1,  1,  0))/sqrt(2),    # 1
        Vec{3,T}(( 1,  1,  0))/sqrt(2),    # 2
        Vec{3,T}(( 1, -1,  0))/sqrt(2),    # 3
        Vec{3,T}(( 1, -1,  0))/sqrt(2),    # 4
        Vec{3,T}(( 1,  0,  1))/sqrt(2),    # 5
        Vec{3,T}(( 1,  0,  1))/sqrt(2),    # 6
        Vec{3,T}(( 1,  0, -1))/sqrt(2),    # 7
        Vec{3,T}(( 1,  0, -1))/sqrt(2),    # 8
        Vec{3,T}(( 0,  1,  1))/sqrt(2),    # 9
        Vec{3,T}(( 0,  1,  1))/sqrt(2),    # 10
        Vec{3,T}(( 0,  1, -1))/sqrt(2),    # 11
        Vec{3,T}(( 0,  1, -1))/sqrt(2),    # 12
        )
    slip_directions = SVector{12,Vec{3,T}}(
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 1
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 2
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 3
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 4
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 5
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 6
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 7
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 8
        Vec{3,T}(( 1,  1, -1))/sqrt(3),    # 9
        Vec{3,T}(( 1, -1,  1))/sqrt(3),    # 10
        Vec{3,T}(( 1,  1,  1))/sqrt(3),    # 11
        Vec{3,T}(( 1, -1, -1))/sqrt(3),    # 12
    )
    return slip_planes, slip_directions
end
