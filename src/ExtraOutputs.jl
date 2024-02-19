update_extras!(::AbstractExtraOutput, args...) = nothing   # If not further specialized, do nothing

# Return the solution (X and dRdXᴹ) to the local problem
# and cache when calculating derivatives
mutable struct DiffOutputHelper{TX, T} <: AbstractExtraOutput
    X::TX
    const dRdX_invᴹ::Matrix{T}
    updated::Bool   # Was the output updated in the timestep?
    const ∂R∂ϵᴹ::Matrix{T}
    const ∂R∂ⁿsᴹ::Matrix{T}
    const ∂R∂pᴹ::Matrix{T}
    const ∂s∂Xᴹ::Matrix{T}
    const ∂s∂ⁿsᴹ::Matrix{T}
    const ∂s∂pᴹ::Matrix{T}
end

function DiffOutputHelper(m::AbstractMaterial)
    T = MMB.get_parameter_type(m)
    X = initial_guess(m, initial_material_state(m), zero(SymmetricTensor{2,3}))
    NR = Tensors.n_components(typeof(X))
    Nσ = 6
    Ns = get_num_statevars(m)
    Np = get_num_params(m)
    dRdX_invᴹ  = zeros(T, NR, NR)
    ∂R∂ϵᴹ  = zeros(T, NR, Nσ)
    ∂R∂ⁿsᴹ = zeros(T, NR, Ns)
    ∂R∂pᴹ  = zeros(T, NR, Np)
    ∂s∂Xᴹ  = zeros(T, Ns, NR)
    ∂s∂ⁿsᴹ = zeros(T, Ns, Ns)
    ∂s∂pᴹ  = zeros(T, Ns, Np)
    return DiffOutputHelper(X, dRdX_invᴹ, false, ∂R∂ϵᴹ, ∂R∂ⁿsᴹ, ∂R∂pᴹ, ∂s∂Xᴹ, ∂s∂ⁿsᴹ, ∂s∂pᴹ)
end

update_extras!(extras::DiffOutputHelper) = (extras.updated = false)

function update_extras!(extras::DiffOutputHelper, X, dRdX_invᴹ)
    extras.X = X
    extras.dRdX_invᴹ .= dRdX_invᴹ
    extras.updated = true
end