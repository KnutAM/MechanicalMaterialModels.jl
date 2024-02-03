function test_numerical_derivative(deriv, deriv_num, m, ϵ, state, Δt, cache, deriv_extra_output, dσdϵ)
    obtain_numerical_material_derivative!(deriv_num, m, ϵ, state, Δt)
    differentiate_material!(deriv, m, ϵ, state, Δt, cache, deriv_extra_output, dσdϵ)
    for field in fieldnames(typeof(deriv))
        @test isapprox(
            getfield(deriv, field),
            getfield(deriv_num, field);
            rtol=1.e-3
            )
    end
end

@testset "Differentiate" begin
    # Create material and vector
    E = 210.e3
    e = LinearElastic(E=E, ν=0.3)
    σ_y0=150.0
    iso = (Voce(Hiso=20.e3, κ∞=200.0), Swift(K=10.0, λ0=1.e-4, n=0.4))
    kin = (ArmstrongFrederick(Hkin=5.e3, β∞=100.0), Delobelle(Hkin=3.e3, β∞=200.0, δ=0.4))
    m1 = Plastic(elastic=e, yield=σ_y0, isotropic=iso, kinematic=kin, overstress=RateIndependent())
    
    for m in (e, m1)
        ## Preallocations
        deriv = MechanicalMaterialModels.MaterialDerivatives(m)
        deriv_num = MechanicalMaterialModels.MaterialDerivatives(m)
        deriv_extra_output = allocate_differentiation_output(m)
        cache = allocate_material_cache(m)
        
        # Test initial derivative
        ## Conditions
        state = initial_material_state(m)
        ϵ = zero(SymmetricTensor{2,3})
        Δt = 1.0
        σ, dσdϵ, new_state = material_response(m, ϵ, state, Δt, cache, deriv_extra_output)
        @testset "Initial ($(nameof(typeof(m))))" begin
            test_numerical_derivative(deriv, deriv_num, m, ϵ, state, Δt, cache, deriv_extra_output, dσdϵ)
        end

        # Test small increase in strain, but still in the elastic regime
        ϵ = SymmetricTensor{2,3}((i,j) -> i==j==1 ? 0.1*σ_y0/E : 0.0)
        σ, dσdϵ, new_state = material_response(m, ϵ, state, Δt, cache, deriv_extra_output)
        @testset "Elastic ($(nameof(typeof(m))))" begin
            test_numerical_derivative(deriv, deriv_num, m, ϵ, state, Δt, cache, deriv_extra_output, dσdϵ)
        end

        # Run until into the plastic regime
        ϵ21 = 0.01
        N = 100
        s21, σ, dσdϵ, state, ϵ = run_shear(m, ϵ21, N)
        ϵ = ϵ*(N+1)/N
        σ, dσdϵ, new_state = material_response(m, ϵ, state, Δt, cache, deriv_extra_output)
        
        @testset "Plastic ($(nameof(typeof(m))))" begin
            test_numerical_derivative(deriv, deriv_num, m, ϵ, state, Δt, cache, deriv_extra_output, dσdϵ)
        end
    end
end

