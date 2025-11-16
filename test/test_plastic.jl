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
    pisov, nisop = MechMat.fromvectortuple(piso, iso)
    @test iso == pisov
    @test nisop == length(piso)
    # Test for the full material
    param = vcat(pel, σ_y0, piso, pkin)
    @test m1 == fromvector(param, m1)
    pvec = zeros(get_vector_length(m1))
    @test length(pvec) == length(param)
    tovector!(pvec, m1)
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

    # Linear hardening
    mlh = Plastic(elastic=e, yield=VonMises(σ_y0), isotropic=Voce(Hiso=Hiso, κ∞=Inf), kinematic=ArmstrongFrederick(Hkin=Hk1, β∞=Inf))
    ϵy = σ_y0/E
    Δϵ = ϵy/N
    s11, σ, dσdϵ, state, ϵ = run_normal(mlh, Δϵ*(N+1), N+1; stress_state=UniaxialStress())
    H = Hk1 + Hiso
    Hp = E*H/(E+H)
    @test s11[end-1] ≈ σ_y0
    @test dσdϵ[1,1,1,1] ≈ Hp 
    @test (s11[end] - s11[end-1])/Δϵ ≈ Hp

    # Non-convergent material, check error throws correctly
    m2 = Plastic(elastic=e, yield=-σ_y0, isotropic=(Voce(Hiso=0.0, κ∞=1.0),), kinematic=(ArmstrongFrederick(Hkin=-10*E, β∞=-1000*E),))
    @test_throws NoLocalConvergence run_shear(m2, Δϵ21, N)

    # Show (note, different conditionals for 1 or more of iso/kin hard)
    show_str = show_as_string(m1)
    @test contains(show_str, "Plastic")
    @test contains(show_str, "LinearElastic")
    @test contains(show_str, "VonMises")
    @test contains(show_str, "Voce")
    @test contains(show_str, "Swift")
    @test contains(show_str, "ArmstrongFrederick")
    @test contains(show_str, "Delobelle")
    @test contains(show_str, "Rate independent")

    show_str = show_as_string(Plastic(;elastic=e, yield=DruckerPrager(1.0, 0.1), isotropic=Voce(0.0, 1.0), kinematic=ArmstrongFrederick(0.0, 1.0), overstress=NortonOverstress(1.0, 1.0)))
    @test contains(show_str, "Plastic")
    @test contains(show_str, "LinearElastic")
    @test contains(show_str, "DruckerPrager")
    @test contains(show_str, "Voce")
    @test contains(show_str, "ArmstrongFrederick")
    @test contains(show_str, "Norton overstress")
end

@testset "SimplePlastic" begin
    # Run an example, compare with Plastic
    E = 210.e3;     ν = 0.3
    Y0 = 150.0
    Hkin = 1 * 10e3;    β∞ = 25.0
    Hiso = 1 * 13e3;    κ∞ = 35.0
    G = E/(2 * (1 + ν))
    K = E / (3 * (1 - 2 * ν))
    elastic = LinearElastic(;E, ν)
    kinematic = ArmstrongFrederick(;Hkin, β∞)
    isotropic = Voce(;Hiso, κ∞)
    m0 = Plastic(;elastic, yield = Y0, kinematic, isotropic)
    m0_simple = SimplePlastic(; G, K, Y0, Hiso, κ∞, Hkin, β∞)
    m0_simple_from_plastic = SimplePlastic(m0)
    for key in fieldnames(typeof(m0_simple))
        @test getproperty(m0_simple, key) ≈ getproperty(m0_simple_from_plastic, key)
    end

    Δϵ21 = 0.01
    N = 100
    s21, σ, dσdϵ, state, ϵ0 = run_shear(m0, Δϵ21, N)

    ϵp_old = state.ϵp
    κold = state.κ[1]
    βold = state.β[1]
    old = MechanicalMaterialModels.SimplePlasticState(ϵp_old, βold, κold)

    ϵt = SymmetricTensor{2,3}((i,j) -> i==2 && j==1 ? ϵ0[i, j] + 1e-4 : 0.0)

    extras = allocate_differentiation_output(m0)
    σ1, dσdϵ1, state1 = material_response(m0, ϵt, state, 0.1, allocate_material_cache(m0), extras)

    σ_s, dσdϵ_s, state_s = material_response(m0_simple, ϵt, old, 0.1)

    Δλ = extras.X.Δλ
    # κ = MechanicalMaterialModels.calculate_κ(Δλ, m0_simple, old)
    σdev, β, κ = MechanicalMaterialModels.calculate_evolution(Δλ, m0_simple, ϵt, old)

    @test state1.κ[1] ≈ state_s.κ ≈ κ
    @test state1.β[1] ≈ state_s.β ≈ β
    @test state1.ϵp ≈ state_s.ϵp
    @test dev(σ1) ≈ dev(σ_s) ≈ σdev
    @test σ1 ≈ σ_s

    Δϵ21 = 0.01
    N = 100
    s21, σ, dσdϵ, state, ϵ = run_shear(m0, Δϵ21, N)
    s21_s, σ_s, dσdϵ_s, state_s, ϵ_s = run_shear(m0_simple, Δϵ21, N)
    # Using tol=1e-10 in newtonsolve inside material_response makes them pass standard ≈ test...
    @test isapprox(s21, s21_s; atol = abs(s21[end]) * 1e-6)
    @test σ ≈ σ_s
    @test dσdϵ ≈ dσdϵ_s
    @test ϵ ≈ ϵ_s

    # Test stability
    γv = [0.0, 0.5, -0.5, 0.0, -1.0, 1.0]
    @test_throws NoLocalConvergence run_shear(m0, γv)
    run_shear(m0_simple, γv)
end