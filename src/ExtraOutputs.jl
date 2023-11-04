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
    const ∂s∂pᴹ::Matrix{T}
end

function DiffOutputHelper(m::Plastic)
    T = get_base_numbertype(m)
    X = initial_guess(m, initial_material_state(m), zero(SymmetricTensor{2,3}))
    NR = Tensors.n_components(typeof(X))
    Nσ = 6
    Ns = get_num_statevars(m)
    Np = get_num_params(m)
    dRdXᴹ  = zeros(T, NR, NR)
    ∂R∂ϵᴹ  = zeros(T, NR, Nσ)
    ∂R∂ⁿsᴹ = zeros(T, NR, Ns)
    ∂R∂pᴹ  = zeros(T, NR, Np)
    ∂s∂Xᴹ  = zeros(T, Ns, NR)
    ∂s∂ⁿsᴹ = zeros(T, Ns, Ns)
    ∂s∂pᴹ  = zeros(T, Ns, Np)
    return DiffOutputHelper(X, dRdXᴹ, false, ∂R∂ϵᴹ, ∂R∂ⁿsᴹ, ∂R∂pᴹ, ∂s∂Xᴹ, ∂s∂ⁿsᴹ, ∂s∂pᴹ)
end

update_extras!(extras::DiffOutputHelper) = (extras.updated = false)

function update_extras!(extras::DiffOutputHelper, X, dRdXᴹ)
    extras.X = X
    extras.dRdXᴹ .= dRdXᴹ
    extras.updated = true
end