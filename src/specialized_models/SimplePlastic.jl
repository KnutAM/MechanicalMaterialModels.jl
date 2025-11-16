@kwdef struct SimplePlastic{T, T_tol} <: AbstractMaterial
    G::T
    K::T 
    Y0::T 
    Hiso::T 
    κ∞::T 
    Hkin::T 
    β∞::T
    maxiter::Int = 100
    tolerance::T_tol = sqrt(eps(typeof(G)))
end

function SimplePlastic(mp::Plastic{Ela, Yld, Iso, Kin, Rat}) where {
        Ela <: LinearElastic{<:Any, :isotropic},
        Yld <: VonMises,
        Iso <: Tuple{<:Voce},
        Kin <: Tuple{<:ArmstrongFrederick},
        Rat <: RateIndependent
        }
    
    E, ν = mp.elastic.p 
    return SimplePlastic(;
        G  = E/(2 * (1 + ν)), K = E / (3 * (1 - 2 * ν)), Y0 = mp.yield.Y0, 
        Hiso = only(mp.isotropic).Hiso, κ∞ = only(mp.isotropic).κ∞,
        Hkin = only(mp.kinematic).Hkin, β∞ = only(mp.kinematic).β∞,
        maxiter = mp.maxiter, tolerance = mp.tolerance)
end

struct SimplePlasticState{T} <: AbstractMaterialState
    ϵp::SymmetricTensor{2,3,T,6}
    β::SymmetricTensor{2,3,T,6}
    κ::T
end

function MMB.initial_material_state(::SimplePlastic)
    return SimplePlasticState(zero(SymmetricTensor{2,3}), zero(SymmetricTensor{2,3}), 0.0)
end

function MMB.material_response(m::SimplePlastic, ϵ::SymmetricTensor{2,3}, old::SimplePlasticState, Δt, args...)
    σ_trial = 2 * m.G * dev(ϵ - old.ϵp) + 3 * m.K * vol(ϵ)
    ϕ_trial = vonmises(σ_trial - old.β) - (m.Y0 + old.κ)
    if ϕ_trial ≤ 0
        I2 = one(ϵ)
        IxI = I2 ⊗ I2
        I4 = one(SymmetricTensor{4,3})
        E4 = 2 * m.G * (I4 - IxI/3) + m.K * IxI
        return σ_trial, E4, old 
    else
        rf(x) = residual(x, m, ϵ, old)
        Δλ, ∂r∂x, converged = newtonsolve(rf, 0.0; maxiter = m.maxiter, tol = m.tolerance)
        #=
        # Using bisection
        Δλ1, ∂r∂x1, converged = newtonsolve(rf, 0.0)
        # Workaround to avoid Δλ and ∂r∂x as `Core.Box`:ed
        Δλ, ∂r∂x = if !converged || Δλ1 < 0 
            Δλ2 = bisect_solve(rf, ϕ_trial / m.G)
            ∂r∂x2 = ForwardDiff.derivative(rf, Δλ2)
            (Δλ2, ∂r∂x2)
        else
            (Δλ1, ∂r∂x1)
        end
        converged = true
        =#

        if converged
            ∂σ∂ϵ = gradient(e -> calculate_sigma(Δλ, m, e, old), ϵ)
            ∂σ∂x, σ = gradient(x -> calculate_sigma(x, m, ϵ, old), Δλ, :all)
            ∂r∂ϵ = gradient(e -> residual(Δλ, m, e, old), ϵ)
            # drdϵ = 0 = ∂r∂ϵ + ∂r∂x * dxdϵ => dxdϵ = - inv(∂r∂x) * ∂r∂ϵ
            dσdϵ = ∂σ∂ϵ - (∂σ∂x / ∂r∂x) ⊗ ∂r∂ϵ

            σdev, β, κ = calculate_evolution(Δλ, m, ϵ, old)
            ϵp = old.ϵp + Δλ * (3/2) * (σdev - β) / (m.Y0 + κ)

            return σ, dσdϵ, SimplePlasticState(ϵp, β, κ)

        else
            throw(MMB.NoLocalConvergence("$(typeof(m)): newtonsolve! didn't converge, ϵ = ", ϵ))
        end
    end
end

function calculate_κ(Δλ, m::SimplePlastic, old::SimplePlasticState)
    return m.κ∞ * (old.κ + Δλ * m.Hiso) / (m.κ∞ + Δλ * m.Hiso)
end

function dβ_fun(Δλ, m::SimplePlastic, κ)
    return m.Y0 + κ + Δλ * m.Hkin * (m.β∞ + m.Y0 + κ) / m.β∞
end

function calculate_sigma(Δλ, m::SimplePlastic, ϵ::SymmetricTensor, old::SimplePlasticState)
    σdev, _, _ = calculate_evolution(Δλ, m, ϵ, old)
    return σdev + 3 * m.K * vol(ϵ)
end

function calculate_evolution(Δλ, m::SimplePlastic, ϵ::SymmetricTensor, old::SimplePlasticState)
    κ = calculate_κ(Δλ, m, old)
    dβ = dβ_fun(Δλ, m, κ)
    Y = m.Y0 + κ
    σdev = (2 * m.G * dev(ϵ) - 2 * m.G * (old.ϵp - (3 * Δλ / (2 * Y)) * (old.β * Y / dβ))) / (
        1 + (3 * m.G * Δλ / Y) * (1 - Δλ * m.Hkin / dβ) )
    β = (old.β * Y + Δλ * m.Hkin * σdev) / dβ
    return σdev, β, κ
end

function residual(Δλ, m::SimplePlastic, ϵ::SymmetricTensor, old::SimplePlasticState)
    σdev, β, κ = calculate_evolution(Δλ, m, ϵ, old)
    Φ = vonmises(σdev - β) - (m.Y0 + κ)
    return Φ
end

#=
# Extra stability, but doesn't seem to be hit in normal cases.
# Therefore, excluded as it is not tested.
function bisect_solve(rf, x0::Number; maxiter = 200, tol = 1e-6)
    x_left = zero(x0)
    x_right = x0
    r_left = rf(x_left)
    @show r_left
    @assert r_left < 0 # No need for premature generalization...

    # Move right untill r_right is above zero 
    r_right = rf(x0)
    iter = 0
    while r_right < 0
        x_left = x_right
        x_right = 2*x_right 
        r_right = rf(x_right)
        println(iter, ": ", (r_right, x_right), " (move right)")
        iter += 1
        iter > maxiter && error("Didn't converge (could not find positive value)")
    end

    # Bisection
    for i in iter:maxiter
        # 0 = r_center ≈ r_left + (r_right - r_left) * (x_center - x_left) / (x_right - x_left)
        x_center = x_left - r_left * (x_right - x_left) / (r_right - r_left)
        r_center = rf(x_center)
        println(iter, ": ", (r_center, x_center), " (bisect)")
        abs(r_center) < tol && return x_center 
        if r_center < 0
            x_left = x_center
            r_left = r_center
        else
            x_right = x_center 
            r_right = r_center 
        end
    end
    error("Did not find sufficiently small residual")
end
=#