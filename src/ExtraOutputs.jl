update_extras!(::AbstractExtraOutput, args...) = nothing   # If not further specialized, do nothing

# Return the solution (X and dRdXᴹ) to the local problem
# and cache when calculating derivatives
mutable struct DiffOutputHelper{TX, T} <: AbstractExtraOutput
    X::TX
    const dRdXᴹ::Matrix{T}
    updated::Bool   # Was the output updated in the timestep?
    const ∂R∂ϵᴹ::Matrix{T}
    const ∂R∂ⁿsᴹ::Matrix{T}
    const ∂R∂pᴹ::Matrix{T}
    const ∂s∂Xᴹ::Matrix{T}
    const ∂s∂ⁿsᴹ::Matrix{T}
end

function DiffOutputHelper(m::Plastic)
    X = initial_guess(m, initial_material_state(m), zero(SymmetricTensor{2,3}))
    NR = Tensors.n_components(typeof(X))
    Nσ = 6
    Ns = get_num_statevars(m)
    Np = get_num_params(m)
    dRdXᴹ = zeros(NR,NR)
    ∂R∂ϵᴹ = zeros(NR,Nσ)
    ∂R∂ⁿsᴹ = zeros(NR,Ns)
    ∂R∂pᴹ = zeros(NR,Np)
    ∂s∂Xᴹ = zeros(Ns,NR)
    ∂s∂ⁿsᴹ = zeros(Ns,Ns)
    return DiffOutputHelper(X, dRdXᴹ, false, ∂R∂ϵᴹ, ∂R∂ⁿsᴹ, ∂R∂pᴹ, ∂s∂Xᴹ, ∂s∂ⁿsᴹ)
end

update_extras!(extras::DiffOutputHelper) = (extras.updated = false)

function update_extras!(extras::DiffOutputHelper, X, dRdXᴹ)
    extras.X = X
    extras.dRdXᴹ .= dRdXᴹ
    extras.updated = true
end