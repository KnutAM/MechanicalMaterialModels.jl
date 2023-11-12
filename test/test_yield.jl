@testset "YieldCriterion" begin
    s = rand()
    σ_normal = SymmetricTensor{2,3}((i,j)->(i==j==1)      ? s : 0.0)
    σ_shear  = SymmetricTensor{2,3}((i,j)->(i==3 && j==1) ? s/√3 : 0.0)
    σ_rand = rand(SymmetricTensor{2,3})
    
    @testset "VonMises" begin
        sy = rand()
        yc = VonMises(sy)
        @test MechMat.effective_stress(yc, σ_normal) ≈ s
        @test MechMat.effective_stress(yc, σ_shear)  ≈ s
        @test MechMat.effective_stress_gradient(yc, σ_rand) ≈ (3/2)*dev(σ_rand)/MechMat.effective_stress(yc, σ_rand)
        # Homogeneous of order 1
        @test MechMat.effective_stress(yc, π*σ_rand) ≈ π*MechMat.effective_stress(yc, σ_rand)

        # Full yield surface 
        k = rand()
        @test MechMat.yield_criterion(yc, σ_shear, k) ≈ (s - (sy+k))

        # Conversions
        @test MMB.get_num_params(yc) == 1
        test_conversion(yc)

        # Show 
        @test contains(show_as_string(yc), "VonMises with")
    end

    @testset "DruckerPrager" begin
        Y0, B = rand(2)
        yc = DruckerPrager(;Y0, B)
        yc_vm = VonMises(Y0)
        # Check same for deviatoric stress
        @test MechMat.effective_stress(yc, dev(σ_rand)) ≈ MechMat.effective_stress(yc_vm, dev(σ_rand))
        # Check direction 
        ν_vm = (3/2)*dev(σ_rand)/MechMat.effective_stress(yc_vm, σ_rand)
        @test MechMat.effective_stress_gradient(yc, σ_rand) ≈ ν_vm - B*one(ν_vm)
        # Homogeneous of order 1
        @test MechMat.effective_stress(yc, π*σ_rand) ≈ π*MechMat.effective_stress(yc, σ_rand)
        
        # Full yield surface 
        k = rand()
        @test MechMat.yield_criterion(yc, σ_shear, k) ≈ (s - (Y0+k))

        # Conversions
        @test MMB.get_num_params(yc) == 2
        test_conversion(yc)

        # Show 
        @test contains(show_as_string(yc), "DruckerPrager with")
    end
end