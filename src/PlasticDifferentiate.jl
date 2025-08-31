MMB.allocate_differentiation_output(m::Plastic) = DiffOutputHelper(m)

function MMB.differentiate_material!(deriv::MaterialDerivatives, m::Plastic, ϵ, ⁿs, Δt, cache, diff_helper::DiffOutputHelper, dσdϵ)
    if diff_helper.updated   # Plastic response
        differentiate_material_plastic!(deriv, m, ϵ, ⁿs, Δt, cache, diff_helper, dσdϵ)
    else                        # Elastic response
        differentiate_material_elastic!(deriv, m, ϵ, ⁿs, diff_helper, dσdϵ)
    end
end

function differentiate_material_elastic!(deriv::MaterialDerivatives{T}, m, ϵ, ⁿs, dh::DiffOutputHelper, dσdϵ) where{T}
    p = tovector(m)
    sv = tovector(ⁿs)
    # Differentiate stress
    tomandel!(deriv.dσdϵ, dσdϵ) # Already calculated

    ## Specialized function to only calculate the elastic stress
    calculate_elastic_stress(m, ϵ, old) = calculate_stress(m.elastic, ϵ - old.ϵp)
    
    ## Differentiate stress by old state   
    σ_from_state(old_vector) = tomandel(calculate_elastic_stress(m, ϵ, fromvector(old_vector, ⁿs)))
    ForwardDiff.jacobian!(dh.dσdⁿsᴹ, σ_from_state, sv)

    ## Differentiate stress by parameters
    σ_from_param(p_vector) = tomandel(calculate_elastic_stress(fromvector(p_vector, m), ϵ, ⁿs))
    ForwardDiff.jacobian!(deriv.dσdp, σ_from_param, p)
    deriv.dσdp .+= dh.dσdⁿsᴹ * deriv.dsdp
    
    # Differentiate state
    fill!(deriv.dsdϵ, zero(T)) # Constant state
    # deriv.dsdp               # Constant state => s depends equally on p as it did in the last time step
end
    
function differentiate_material_plastic!(deriv::MaterialDerivatives, m::Plastic, ϵ, ⁿs, Δt, cache, diff_helper::DiffOutputHelper, dσdϵ)
    # Extract from input
    p = tovector(m)
    sv = tovector(ⁿs)
    X = diff_helper.X
    ∂R∂Xinvᴹ = diff_helper.dRdX_invᴹ
    
    Nσ = size(deriv.dσdϵ, 1)
    ∂R∂ϵᴹ = diff_helper.∂R∂ϵᴹ
    ∂R∂ⁿsᴹ = diff_helper.∂R∂ⁿsᴹ
    ∂R∂pᴹ = diff_helper.∂R∂pᴹ
    ∂s∂Xᴹ = diff_helper.∂s∂Xᴹ
    ∂s∂pᴹ = diff_helper.∂s∂pᴹ
    dσdⁿsᴹ = diff_helper.dσdⁿsᴹ
    dsdⁿsᴹ = diff_helper.dsdⁿsᴹ
    
    # Precalculations
    R_from_strain(ϵ_vector) = tovector(residual(X, m, ⁿs, frommandel(baseof(ϵ), ϵ_vector), Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂ϵᴹ, R_from_strain, tomandel(ϵ))

    R_from_state(old_vector) = tovector(residual(X, m, fromvector(old_vector, ⁿs), ϵ, Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂ⁿsᴹ, R_from_state, sv)

    R_from_param(p_vector) = tovector(residual(X, fromvector(p_vector, m), ⁿs, ϵ, Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂pᴹ, R_from_param, p)

    # Differentiate the stress
    ## dσdXᴹ = [I 0] where I is 6x6 identity matrix
    ∂σ∂Xᴹ_times_∂R∂Xinvᴹ = ∂R∂Xinvᴹ[1:Nσ, :]
    ## dσdϵ
    tomandel!(deriv.dσdϵ, dσdϵ) # Already calculated
    ## dσdⁿs
    dσdⁿsᴹ .= -∂σ∂Xᴹ_times_∂R∂Xinvᴹ * ∂R∂ⁿsᴹ # ∂σ∂ⁿs=0 (σ=X.σ)
    ## dσdp
    deriv.dσdp .= -∂σ∂Xᴹ_times_∂R∂Xinvᴹ * ∂R∂pᴹ # ∂σ∂p=0 (σ=X.σ)
    deriv.dσdp .+= dσdⁿsᴹ * deriv.dsdp

    # Differentiate the new state 
    ## Calculate ∂s∂Xᴹ first
    s_from_X(X_vector) = tovector(get_plastic_result(m, fromvector(X_vector, X), ⁿs)[2])
    ForwardDiff.jacobian!(∂s∂Xᴹ, s_from_X, tovector(X))
    ∂s∂Xᴹ_times_∂R∂Xinvᴹ = ∂s∂Xᴹ*∂R∂Xinvᴹ

    ## dsdϵ
    deriv.dsdϵ .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂ϵᴹ    # ∂s∂ϵ = 0 (get_plastic_result does not take the strain as argument)

    ## dsdⁿs
    s_from_state(old_vector) = tovector(get_plastic_result(m, X, fromvector(old_vector, ⁿs))[2])
    ForwardDiff.jacobian!(dsdⁿsᴹ, s_from_state, sv)
    dsdⁿsᴹ .-= ∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂ⁿsᴹ

    ## dsdp
    s_from_p(p_vector) = tovector(get_plastic_result(fromvector(p_vector, m), X, ⁿs)[2])
    ForwardDiff.jacobian!(∂s∂pᴹ, s_from_p, p)
    deriv.dsdp .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂pᴹ .+ ∂s∂pᴹ .+ dsdⁿsᴹ * deriv.dsdp
end
