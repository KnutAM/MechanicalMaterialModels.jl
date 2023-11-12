@testset "YieldCriterion" begin
    s = rand()
    σ_normal = SymmetricTensor{2,3}((i,j)->(i==j==1)      ? s : 0.0)
    σ_shear  = SymmetricTensor{2,3}((i,j)->(i==3 && j==1) ? s : 0.0)
    σ_rand = rand(SymmetricTensor{2,3})
    
    @testset "VonMises" begin
        sy = rand()
        yc = VonMises(sy)
        @test MechMat.effective_stress(yc, σ_normal) ≈ s
        @test MechMat.effective_stress(yc, σ_shear) ≈ √3*s
        @test MechMat.effective_stress_gradient(yc, σ_rand) ≈ (3/2)*dev(σ_rand)/MechMat.effective_stress(yc, σ_rand)
        # Homogeneous of order 1
        @test MechMat.effective_stress(yc, π*σ_rand) ≈ π*MechMat.effective_stress(yc, σ_rand)
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
    end
end