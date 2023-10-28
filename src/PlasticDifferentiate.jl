MMB.allocate_differentiation_output(m::Plastic) = DiffOutputHelper(m)
function MMB.differentiate_material!(deriv::MaterialDerivatives, m::Plastic, ϵ, ⁿs, Δt, dσdϵ, diff_helper::DiffOutputHelper, cache)
    if diff_helper.updated   # Plastic response
        differentiate_material_plastic!(deriv, m, ϵ, ⁿs, Δt, dσdϵ, diff_helper, cache)
    else                        # Elastic response
        differentiate_material_elastic!(deriv, m, ϵ, ⁿs, dσdϵ)
    end
end

function differentiate_material_elastic!(deriv::MaterialDerivatives{T}, m::Plastic, ϵ, ⁿs, dσdϵ) where{T}
    p = material2vector(m)
    sv = zeros(Tensors.n_components(typeof(ⁿs)))
    # Differentiate stress
    tomandel!(deriv.dσdϵ, dσdϵ) # Already calculated

    ## Specialized function to only calculate the elastic stress
    calculate_elastic_stress(m::Plastic, ϵ, old::PlasticState) = calculate_stress(m.elastic, ϵ-old.ϵₚ)
    
    ## Differentiate stress by old state   
    σ_from_state(old_vector) = tomandel(calculate_elastic_stress(m, ϵ, frommandel(typeof(ⁿs), old_vector)))
    ForwardDiff.jacobian!(deriv.dσdⁿs,σ_from_state,tomandel!(sv,ⁿs))

    ## Differentiate stress by parameters
    σ_from_param(p_vector) = tomandel(calculate_elastic_stress(vector2material(p_vector, m), ϵ, ⁿs))
    ForwardDiff.jacobian!(deriv.dσdp,σ_from_param,p)
    
    # Differentiate state
    fill!(deriv.dsdϵ,  zero(T)) # Constant state
    # new state equal to old state => derivative is the identity matrix
    copyto!(deriv.dsdⁿs, LinearAlgebra.I)
    fill!(deriv.dsdp, zero(T)) # Constant state
end
    
function differentiate_material_plastic!(deriv::MaterialDerivatives, m::Plastic, ϵ, ⁿs, Δt, dσdϵ, diff_helper::DiffOutputHelper, cache)
    # Extract from input
    p = material2vector(m)
    X = diff_helper.X
    ∂R∂Xᴹ = diff_helper.dRdXᴹ
    
    # Preallocations (could be done beforehand also) (X=Nx, P=Np)
    Nσ = size(deriv.dσdϵ, 1)
    ∂R∂ϵᴹ = diff_helper.∂R∂ϵᴹ
    ∂R∂ⁿsᴹ = diff_helper.∂R∂ⁿsᴹ
    ∂R∂pᴹ = diff_helper.∂R∂pᴹ
    ∂s∂Xᴹ = diff_helper.∂s∂Xᴹ
    ∂s∂ⁿsᴹ = diff_helper.∂s∂ⁿsᴹ
    
    # Precalculations
    ⁿs_vector = tomandel(ⁿs)
    ∂R∂Xinvᴹ = inv(∂R∂Xᴹ)

    R_from_strain(ϵ_vector) = tomandel(residual(X, m, ⁿs, frommandel(shapeof(ϵ), ϵ_vector), Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂ϵᴹ, R_from_strain, tomandel(ϵ))

    R_from_state(old_vector) = tomandel(residual(X, m, frommandel(typeof(ⁿs), old_vector), ϵ, Δt, cache.resid))
    ForwardDiff.jacobian!(∂R∂ⁿsᴹ, R_from_state, ⁿs_vector)

    R_from_param(p_vector) = tomandel(residual(X, vector2material(p_vector, m), ⁿs, ϵ, Δt, cache.resid))
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
    s_from_X(X_vector) = tomandel(get_plastic_result(frommandel(typeof(X), X_vector), ⁿs)[2])
    ForwardDiff.jacobian!(∂s∂Xᴹ, s_from_X, tomandel(X))
    ∂s∂Xᴹ_times_∂R∂Xinvᴹ = ∂s∂Xᴹ*∂R∂Xinvᴹ
    ## dsdϵ
    deriv.dsdϵ .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂ϵᴹ    # ∂s∂ϵ = 0 (get_plastic_result does not take the strain as argument)
    ## dsdⁿs
    s_from_state(old_vector) = tomandel(get_plastic_result(X, frommandel(typeof(ⁿs), old_vector))[2])
    ForwardDiff.jacobian!(∂s∂ⁿsᴹ, s_from_state, ⁿs_vector)
    deriv.dsdⁿs .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂ⁿsᴹ + ∂s∂ⁿsᴹ
    ## dsdp
    deriv.dsdp .= -∂s∂Xᴹ_times_∂R∂Xinvᴹ*∂R∂pᴹ    # ∂s∂p = 0 (get_plastic_result does not take the material as argument)

end