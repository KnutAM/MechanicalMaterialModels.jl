function test_numerical_derivative(deriv, deriv_num, m, ϵ, state, Δt, cache, deriv_extra_output, dσdϵ)
    MatTest.obtain_numerical_material_derivative!(deriv_num, m, ϵ, state, Δt)
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

@testset "Numerical (3D material)" begin
        E = 210.e3
        ν = 0.3
        e = LinearElastic(;E, ν)
        σ_y0=150.0
        iso = (Voce(Hiso=20.e3, κ∞=200.0), Swift(K=10.0, λ0=1.e-4, n=0.4))
        kin = (ArmstrongFrederick(Hkin=5.e3, β∞=100.0), Delobelle(Hkin=3.e3, β∞=200.0, δ=0.4))
        m1 = Plastic(elastic=e, yield=σ_y0, isotropic=iso, kinematic=kin, overstress=RateIndependent())
        m2 = Plastic(elastic=e, yield=σ_y0, isotropic=iso, kinematic=kin, overstress=NortonOverstress(; tstar = 0.123, nexp = 2.5))
        materials = ["elastic" => e, "plastic" => m1, "viscoplastic" => m2]
        for (key, m) in materials
            @testset "$key" begin
                @testset "Initial response" begin
                    ϵ = rand(SymmetricTensor{2,3}) * 0.1 * σ_y0 / E
                    Δt = 1e-2
                    MatTest.test_derivative(m, ϵ, initial_material_state(m), Δt; 
                        numdiffsettings = (fdtype = Val{:central},),
                        comparesettings = (atol_min = 1e-10, rtol_min = 1e-8),
                        )
                    diff = MaterialDerivatives(m)
                    copy!(diff.dsdp, rand(size(diff.dsdp)...))
                    MatTest.test_derivative(m, ϵ, initial_material_state(m), Δt; 
                        numdiffsettings = (fdtype = Val{:central},),
                        comparesettings = (atol_min = 1e-10, rtol_min = 1e-8),
                        diff)
                end
                @testset "Just past yielding" begin
                    ϵ = SymmetricTensor{2,3}((i,j) -> i==j==1 ? 1.0 : (i == j ? -ν : 0.0)) * 1.01 * σ_y0 / E
                    Δt = 1e-2
                    MatTest.test_derivative(m, ϵ, initial_material_state(m), Δt; 
                        numdiffsettings = (fdtype = Val{:central}, relstep = 1e-8),
                        comparesettings = (atol_min = 1e-8, rtol_min = 1e-4, print_tol = false)
                        )
                end
                @testset "After shear loading" begin
                    num_steps = 10; t_end = 0.01
                    ϵ21 = num_steps * 0.1 * 3 * σ_y0 / E;
                    relstep = 1e-6
                    stressfun(p) = MatTest.runstrain(fromvector(p, m), ϵ21, (2, 1), t_end, num_steps)[1]
                    dσ21_dp_num = FiniteDiff.finite_difference_jacobian(stressfun, BigFloat.(tovector(m)), Val{:central}; relstep)
                    σv, state, dσ21_dp, diff = MatTest.runstrain_diff(m, ϵ21, (2, 1), t_end, num_steps)
                    @test σv ≈ stressfun(tovector(m))
                    #@test dσ21_dp ≈ dσ21_dp_num
                    isapprox(dσ21_dp_num, dσ21_dp; atol = 1000 * eps())
                    scaled_error, maxtol = MatTest.compare_derivatives(dσ21_dp, dσ21_dp_num, σv, tovector(m) * relstep; atol_min = 1e-12, rtol_min = 1e-4, print_tol = false)
                    @test all(x -> x ≤ 1, scaled_error)
                end
                @testset "FullStressState" begin
                    ϵ21 = 3 * σ_y0 / E; num_steps = 10; t_end = 0.01
                    ss = FullStressState()
                    # Check that we get the same result for MatTest.runstresstate and MatTest.runstrain
                    σ_ss, s_ss = MatTest.runstresstate(ss, m, ϵ21, (2, 1), t_end, num_steps)
                    σ, s = MatTest.runstrain(m, ϵ21, (2, 1), t_end, num_steps)
                    @test σ_ss ≈ σ
                    @test tovector(s_ss) ≈ tovector(s)
                end
                for (stress_state, ij) in (
                        (UniaxialStress(), (1,1)), 
                        (UniaxialNormalStress(), (1,2)),
                        (GeneralStressState(SymmetricTensor{2,3,Bool}((true, false, false, false, true, true)), rand(SymmetricTensor{2,3})), (2,2))
                        )
                    @testset "$(nameof(typeof(stress_state))), (i,j) = ($(ij[1]), $(ij[2]))" begin
                        num_steps = 10; t_end = 0.01
                        ϵij = num_steps * 0.1 * 3 * σ_y0 / E
                        relstep = 1e-6
                        stressfun(p) = MatTest.runstresstate(stress_state, fromvector(p, m), ϵij, ij, t_end, num_steps)[1]
                        dσij_dp_num = FiniteDiff.finite_difference_jacobian(stressfun, BigFloat.(tovector(m)), Val{:central}; relstep)
                        σv, state, dσij_dp, diff = MatTest.runstresstate_diff(stress_state, m, ϵij, ij, t_end, num_steps)
                        isapprox(dσij_dp_num, dσij_dp; atol = 1000 * eps())
                        
                        scaled_error, maxtol = MatTest.compare_derivatives(dσij_dp, dσij_dp_num, σv, tovector(m) * relstep; atol_min = 1e-6, rtol_min = 1e-4, print_tol = false)
                        if !all(x -> x ≤ 1, scaled_error)
                            display(scaled_error)
                        end
                        @test all(x -> x ≤ 1, scaled_error)
                    end
                end
            end
        end
    end
    