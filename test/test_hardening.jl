
@testset "KinematicHardening" begin
    get_evolution(laws::Tuple, ν, β) = map(law->MechMat.get_evolution(law, ν, β), laws)

    Hkin=10.e3; β∞=30.0; δ=0.5; m=3.9
    af = ArmstrongFrederick(Hkin=Hkin, β∞=β∞);
    db = Delobelle(Hkin=Hkin, β∞=β∞, δ=δ)
    ow = OhnoWang(Hkin=Hkin, β∞=β∞, m=m)
    σ_red_dev = dev(rand(SymmetricTensor{2,3}))
    ν = (3/2)*σ_red_dev/MechMat.vonmises(σ_red_dev)
    β = zero(SymmetricTensor{2,3})

    # Check that all give the same initial hardening modulus
    af_h0, db_h0, ow_h0 = get_evolution((af, db, ow), ν, β)
    @test af_h0 ≈ db_h0
    @test af_h0 ≈ ow_h0

    # Check that ArmstrongFrederick and Delobelle give the same result when β and ν are aligned
    β = (1.0 + rand())*ν    # Ensure scaling > 0
    af_h1, db_h1, ow_h1 = get_evolution((af, db, ow), ν, β)
    @test af_h1 ≈ db_h1     # Should be equal
    @test !isapprox(af_h0, ow_h1)  # Should not be equal
    
    # Check constructors from vectors
    @test af == fromvector([Hkin, β∞], af)
    @test af == fromvector([1.0, Hkin, β∞], af, offset=1)
    @test db == fromvector([Hkin, β∞, δ], db)
    @test ow == fromvector([Hkin, β∞, m], ow)

    # Test conversions
    for T in (Float64, MatTest.DualT{Float32})
        for hardlaw in (af, db, ow)
            MatTest.test_vectorconversion(T, hardlaw)
        end
    end

    # Check show methods 
    @test contains(show_as_string(af), "ArmstrongFrederick with")
    @test contains(show_as_string(db), "Delobelle with")
    @test contains(show_as_string(ow), "OhnoWang with")
end

@testset "IsotropicHardening" begin
    λ = 0.1
    
    # Voce hardening
    (Hiso, κ∞) = vv = [200.0, 10.0]
    voce = Voce(Hiso=Hiso, κ∞=κ∞)

    # Test conversion from vector
    @test voce == fromvector(vv, voce)
    
    # Calculate Voce response using analytical function
    vocelaw(param, λ) = param.κ∞*(one(λ) - exp(-param.Hiso*λ/param.κ∞))
    κ = vocelaw(voce, λ)
    dκdλ = ForwardDiff.derivative(λarg->vocelaw(voce, λarg), λ)

    @test dκdλ ≈ MechMat.get_evolution(voce, κ)

    # Swift hardening
    K, λ0, n = sv = [10.0, 1.e-3, 2.0]
    swift = Swift(K=K, λ0=λ0, n=n)
    swiftlaw(param, λ) = param.K * (param.λ0 + λ)^param.n 

    # Test conversion from vector
    @test swift == fromvector(sv, swift)

    κ0 = swiftlaw(swift, 0.0)
    κ = swiftlaw(swift, λ)
    dκdλ = ForwardDiff.derivative(λarg->swiftlaw(swift, λarg), λ)
    
    @test dκdλ ≈ MechMat.get_evolution(swift, κ)
    @test κ0 ≈ MechMat.get_initial_value(swift)

        # Test conversions
    for T in (Float64, MatTest.DualT{Float32})
        for hardlaw in (voce, swift)
            MatTest.test_vectorconversion(T, hardlaw)
        end
    end

    # Show
    @test contains(show_as_string(voce), "Voce with")
    @test contains(show_as_string(swift), "Swift with")
end 