@testset "YieldCriterion" begin
    for TB in (SymmetricTensor{2,3}, Tensor{2,3})
        s = rand()
        σ_normal = TB((i,j)->(i==j==1)      ? s : 0.0)
        σ_shear  = TB((i,j)->(i==3 && j==1) || (i==1 && j==3) ? s/√3 : 0.0)
        σ_rand = rand(TB)
        
        @testset "VonMises ($TB)" begin
            sy = rand()
            yc = VonMises(sy)
            @test MechMat.effective_stress(yc, σ_normal) ≈ s
            @test MechMat.effective_stress(yc, σ_shear)  ≈ s
            @test MechMat.effective_stress_gradient(yc, σ_rand) ≈ (3/2)*dev(σ_rand')/MechMat.effective_stress(yc, σ_rand)
            @test MechMat.effective_stress_gradient(yc, σ_rand) ≈ gradient(s->MechMat.effective_stress(yc, s), σ_rand)
            # Homogeneous of order 1
            @test MechMat.effective_stress(yc, π*σ_rand) ≈ π*MechMat.effective_stress(yc, σ_rand)

            # Full yield surface 
            k = rand()
            @test MechMat.yield_criterion(yc, σ_shear, k) ≈ (s - (sy+k))

            # Conversions
            @test get_vector_length(yc) == 1
            MatTest.test_vectorconversion(Float64, yc)
            MatTest.test_vectorconversion(Float32, yc)

            # Show 
            @test contains(show_as_string(yc), "VonMises with")
        end
        
        @testset "DruckerPrager ($TB)" begin
            Y0, B = rand(2)
            yc = DruckerPrager(;Y0, B)
            yc_vm = VonMises(Y0)
            # Check same for deviatoric stress
            @test MechMat.effective_stress(yc, dev(σ_rand)) ≈ MechMat.effective_stress(yc_vm, dev(σ_rand))
            # Check direction (both using AD and analytical, in case analytical is implemented later)
            
            ν_vm = MechMat.effective_stress_gradient(yc_vm, σ_rand)
            @test MechMat.effective_stress_gradient(yc, σ_rand) ≈ ν_vm - B*one(ν_vm)
            @test MechMat.effective_stress_gradient(yc, σ_rand) ≈ gradient(s->MechMat.effective_stress(yc, s), σ_rand)
            # Homogeneous of order 1
            @test MechMat.effective_stress(yc, π*σ_rand) ≈ π*MechMat.effective_stress(yc, σ_rand)
            
            # Full yield surface 
            k = rand()
            @test MechMat.yield_criterion(yc, σ_shear, k) ≈ (s - (Y0+k))

            # Conversions
            @test get_vector_length(yc) == 2
            MatTest.test_vectorconversion(Float64, yc)
            MatTest.test_vectorconversion(MatTest.DualT{Float32}, yc)

            # Show
            @test contains(show_as_string(yc), "DruckerPrager with")
        end
    end

    @testset "DruckerPrager apex" begin
        # At the hydrostatic apex, dev(σred) = 0, which is non-differentiable for the
        # plain vonmises(σred) = sqrt(3/2*dev(σred)⊡dev(σred)) expression (0/0). This must
        # work both for `effective_stress_gradient` itself, and for `effective_stress`
        # when differentiated through (e.g. by a Newton solver computing its own
        # Jacobian via automatic differentiation of the whole residual, which goes
        # through `effective_stress`, not `effective_stress_gradient`).
        for TB in (SymmetricTensor{2,3}, Tensor{2,3}), p in (-1.0, 0.0, 4.0) # exactly representable, so dev(p*I) is exactly 0
            Y0, B = rand(2)
            yc = DruckerPrager(;Y0, B)
            σ_apex = p * one(TB)
            @test iszero(dev(σ_apex))

            ν = MechMat.effective_stress_gradient(yc, σ_apex)
            @test !any(isnan, ν)
            @test ν ≈ -B * one(σ_apex)

            # Differentiating `effective_stress` itself (as a Newton solver's automatic
            # differentiation of its residual would) must not produce NaN/Inf, and must
            # agree with the hand-coded (sub)gradient above.
            ν_ad = gradient(s -> MechMat.effective_stress(yc, s), σ_apex)
            @test !any(isnan, ν_ad)
            @test ν_ad ≈ ν
        end
    end
end