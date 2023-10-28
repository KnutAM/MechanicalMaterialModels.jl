@testset "Plastic" begin
    (E, ν) = pel = [210.e3, 0.3]
    e = LinearElastic(E=E, ν=ν)
    σ_y0=150.0
    (Hiso, κ∞, K, λ0, n) = piso = [20.e3, 200.0, 10.0, 1.e-4, 0.4]
    iso = (Voce(Hiso=Hiso, κ∞=κ∞), Swift(K=K, λ0=λ0, n=n))
    (Hk1, β1, Hk2, β2, δ2) = pkin = [5.e3, 100.0, 3.e3, 200., 0.4]
    kin = (ArmstrongFrederick(Hkin=Hk1, β∞=β1), Delobelle(Hkin=Hk2, β∞=β2, δ=δ2))
    m1 = Plastic(elastic=e, yield=σ_y0, isotropic=iso, kinematic=kin)

    # Test constructors from vectors
    # Single test for just the tuple part
    pisov, nisop = MechMat.vector2materialtuple(piso, iso)
    @test iso == pisov
    @test nisop == length(piso)
    # Test for the full material
    param = vcat(pel, σ_y0, piso, pkin)
    @test m1 == MMB.vector2material(param, m1)
    pvec = zeros(get_num_params(m1))
    @test length(pvec) == length(param)
    MMB.material2vector!(pvec, m1)
    @test pvec == param
    
    # Run an example, not test of result, just checking that it runs
    Δϵ21 = 0.01
    N = 100
    s21, σ, dσdϵ, state, ϵ = run_shear(m1, Δϵ21, N)

    # Ideally plastic, check von-mises stress
    m2 = Plastic(elastic=e, yield=σ_y0, isotropic=(Voce(Hiso=0.0, κ∞=1.0),), kinematic=(ArmstrongFrederick(Hkin=0.0, β∞=1.0),))
    s21, σ, dσdϵ, state, ϵ = run_shear(m2, Δϵ21, N)
    @test MechMat.vonmises(σ) ≈ σ_y0
    s11, σ, dσdϵ, state, ϵ = run_normal(m2, Δϵ21, N)
    @test MechMat.vonmises(σ) ≈ σ_y0

    # Non-convergent material, check error throws correctly
    m2 = Plastic(elastic=e, yield=-σ_y0, isotropic=(Voce(Hiso=0.0, κ∞=1.0),), kinematic=(ArmstrongFrederick(Hkin=-10*E, β∞=-1000*E),))
    @test_throws NoLocalConvergence run_shear(m2, Δϵ21, N)
end
