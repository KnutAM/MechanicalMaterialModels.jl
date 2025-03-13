MMB.allocate_differentiation_output(m::Plastic) = DiffOutputHelper(m)
function MMB.differentiate_material!(deriv::MaterialDerivatives, m::Plastic, ϵ, ⁿs, Δt, cache, diff_helper::DiffOutputHelper, dσdϵ)
    if diff_helper.updated   # Plastic response
        differentiate_material_plastic!(deriv, m, ϵ, ⁿs, Δt, cache, diff_helper, dσdϵ)
    else                        # Elastic response
        differentiate_material_elastic!(deriv, m, ϵ, ⁿs, dσdϵ)
    end
end

function differentiate_material_elastic!(deriv::MaterialDerivatives{T}, m, ϵ, ⁿs, dσdϵ) where{T}
    p = tovector(m)
    sv = zeros(get_num_statevars(ⁿs))
    # Differentiate stress
    tomandel!(deriv.dσdϵ, dσdϵ) # Already calculated

    ## Specialized function to only calculate the elastic stress
    calculate_elastic_stress(m, ϵ, old) = calculate_stress(m.elastic, ϵ - old.ϵp)
    
    ## Differentiate stress by old state   
    σ_from_state(old_vector) = tomandel(calculate_elastic_stress(m, ϵ, fromvector(old_vector, ⁿs)))
    ForwardDiff.jacobian!(deriv.dσdⁿs, σ_from_state, tovector!(sv, ⁿs))

    ## Differentiate stress by parameters
    σ_from_param(p_vector) = tomandel(calculate_elastic_stress(fromvector(p_vector, m), ϵ, ⁿs))
    ForwardDiff.jacobian!(deriv.dσdp, σ_from_param, p)
    
    # Differentiate state
    fill!(deriv.dsdϵ,  zero(T)) # Constant state
    # new state equal to old state => derivative is the identity matrix
    copyto!(deriv.dsdⁿs, LinearAlgebra.I)
    fill!(deriv.dsdp, zero(T)) # Constant state
end
    
function differentiate_material_plastic!(deriv::MaterialDerivatives, m::Plastic, ϵ, ⁿs, Δt, cache, diff_helper::DiffOutputHelper, dσdϵ)
    # Extract from input
    p = tovector(m)
    X = diff_helper.X
    ∂R∂Xinvᴹ = diff_helper.dRdX_invᴹ
    
    Nσ = size(deriv.dσdϵ, 1)
    ∂R∂ϵᴹ = diff_helper.∂R∂ϵᴹ
    ∂R∂ⁿsᴹ = diff_helper.∂R∂ⁿsᴹ
    ∂R∂pᴹ = diff_helper.∂R∂pᴹ
    ∂s∂Xᴹ = diff_helper.∂s∂Xᴹ
    ∂s∂ⁿsᴹ = diff_helper.∂s∂ⁿsᴹ
    ∂s∂pᴹ = diff_helper.∂s∂pᴹ
    
    # Precalculations
    ⁿs_vector = tovector(ⁿs)

    R_from_strain(ϵ_vector) = tovector(residual(X, m, ⁿs, frommandel(baseof(ϵ), ϵ_vector), Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂ϵᴹ, R_from_strain, tomandel(ϵ))

    R_from_state(old_vector) = tovector(residual(X, m, fromvector(old_vector, ⁿs), ϵ, Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂ⁿsᴹ, R_from_state, ⁿs_vector)

    R_from_param(p_vector) = tovector(residual(X, fromvector(p_vector, m), ⁿs, ϵ, Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂pᴹ, R_from_param, p)

    # Differentiate the stress
    ## dσdXᴹ = [I 0] where I is 6x6 identity matrix
    ∂σ∂Xᴹ_times_∂R∂Xinvᴹ = ∂R∂Xinvᴹ[1:Nσ, :]
    ## dσdϵ 
    tomandel!(deriv.dσdϵ, dσdϵ) # Already calculated
    ## dσdⁿs 
    deriv.dσdⁿs .= -∂σ∂Xᴹ_times_∂R∂Xinvᴹ * ∂R∂ⁿsᴹ # ∂σ∂ⁿs=0 (σ=X.σ)
    ## dσdp
    deriv.dσdp .= -∂σ∂Xᴹ_times_∂R∂Xinvᴹ * ∂R∂pᴹ # ∂σ∂p=0 (σ=X.σ)

    # Differentiate the new state 
    ## Calculate ∂s∂Xᴹ first
    s_from_X(X_vector) = tovector(get_plastic_result(m, fromvector(X_vector, X), ⁿs)[2])
    ForwardDiff.jacobian!(∂s∂Xᴹ, s_from_X, tovector(X))
    ∂s∂Xᴹ_times_∂R∂Xinvᴹ = ∂s∂Xᴹ*∂R∂Xinvᴹ
    ## dsdϵ
    deriv.dsdϵ .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂ϵᴹ    # ∂s∂ϵ = 0 (get_plastic_result does not take the strain as argument)
    ## dsdⁿs
    s_from_state(old_vector) = tovector(get_plastic_result(m, X, fromvector(old_vector, ⁿs))[2])
    ForwardDiff.jacobian!(∂s∂ⁿsᴹ, s_from_state, ⁿs_vector)
    deriv.dsdⁿs .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂ⁿsᴹ + ∂s∂ⁿsᴹ
    ## dsdp
    s_from_p(p_vector) = tovector(get_plastic_result(fromvector(p_vector, m), X, ⁿs)[2])
    ForwardDiff.jacobian!(∂s∂pᴹ, s_from_p, p)
    deriv.dsdp .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂pᴹ + ∂s∂pᴹ
    
end