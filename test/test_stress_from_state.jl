@testset "stress_from_state" begin
    @testset "LinearElastic" begin
        m = LinearElastic(E=210.e3, ν=0.3)
        state = initial_material_state(m)
        ϵ = rand(SymmetricTensor{2,3})
        @test stress_from_state(m, ϵ, state) ≈ m.C ⊡ ϵ

        # Reduced stress state: matches material_response's own (iterative) result,
        # and the analytical plane-stress relation
        E, ν = 210.e3, 0.3
        me = LinearElastic(; E, ν)
        rss = ReducedStressState(PlaneStress(), me)
        ϵ11 = 0.01
        ϵ_red = SymmetricTensor{2,2}((ϵ11, 0.0, 0.0))
        σ_direct = stress_from_state(rss, ϵ_red, initial_material_state(rss))
        σ_mr, _, _, _ = material_response(PlaneStress(), me, ϵ_red, initial_material_state(me))
        @test σ_direct ≈ σ_mr
        @test σ_direct[1, 1] ≈ E / (1 - ν^2) * ϵ11
        @test σ_direct[2, 2] ≈ E / (1 - ν^2) * ν * ϵ11
    end

    @testset "Plastic" begin
        E, ν = 210.e3, 0.3
        e = LinearElastic(; E, ν)
        m = Plastic(elastic=e, yield=100.0, isotropic=Voce(Hiso=10.e3, κ∞=200.0), kinematic=ArmstrongFrederick(Hkin=1.e4, β∞=150.0))

        # Load to a converged, plastically loaded state
        state0 = initial_material_state(m)
        ϵ1 = SymmetricTensor{2,3}((i, j) -> (i, j) == (1, 1) ? 0.01 : 0.0)
        σ1, _, state1 = material_response(m, ϵ1, state0, nothing)
        @test stress_from_state(m, ϵ1, state1) ≈ σ1

        # Frozen-state postprocessing: a different strain should give a purely
        # elastic increment from state1, NOT a fresh plastic correction.
        ϵ2 = ϵ1 + SymmetricTensor{2,3}((i, j) -> (i, j) == (1, 1) ? 0.02 : 0.0)
        σ2_frozen = stress_from_state(m, ϵ2, state1)
        @test σ2_frozen ≈ σ1 + e.C ⊡ (ϵ2 - ϵ1)
        σ2_true, _, state2_true = material_response(m, ϵ2, state1, nothing)
        @test !(σ2_true ≈ σ2_frozen) # material_response would further evolve plastically
        @test state2_true.ϵp != state1.ϵp

        # Reduced stress state (mirrors the FerriteAssembly#94 plane-stress postprocessing fix)
        rss = ReducedStressState(PlaneStress(), m)
        ϵ1_red = SymmetricTensor{2,2}((0.01, 0.0, 0.0))
        state0_red = initial_material_state(rss)
        σ1_red, _, state1_red, _ = material_response(rss, ϵ1_red, state0_red, nothing)
        σ1_red_current = stress_from_state(rss, ϵ1_red, state1_red)
        @test σ1_red_current ≈ σ1_red

        # PlaneStrain: verify the transverse plastic strain's elastic coupling is
        # retained (regression check for an issue found in review of an earlier,
        # since-removed Plastic-specific reduced-dimensional method; kept to
        # confirm the generic fallback gets this right too).
        rss_strain = ReducedStressState(PlaneStrain(), m)
        state0_strain = initial_material_state(rss_strain)
        σ1_strain, _, state1_strain, _ = material_response(rss_strain, ϵ1_red, state0_strain, nothing)
        @test state1_strain.ϵp[3, 3] != 0 # sanity: this test only matters if ϵp33 != 0
        σ1_strain_current = stress_from_state(rss_strain, ϵ1_red, state1_strain)
        @test σ1_strain_current ≈ σ1_strain
    end

    @testset "GeneralizedMaxwell" begin
        me = LinearElastic(E=210.e3, ν=0.3)
        chain = Maxwell(G=1.e3, t=1.0)
        m = GeneralizedMaxwell(me, chain)

        state0 = initial_material_state(m)
        ϵ1 = rand(SymmetricTensor{2,3}) / 100
        σ1, _, state1 = material_response(m, ϵ1, state0, 0.5)
        @test stress_from_state(m, ϵ1, state1) ≈ σ1

        # Frozen-state: evaluating at a different strain must not re-solve the
        # viscous strain evolution (which requires Δt); it should be a pure
        # elastic-type increment using the given (fixed) viscous strain.
        ϵ2 = ϵ1 + rand(SymmetricTensor{2,3}) / 100
        σ2_frozen = stress_from_state(m, ϵ2, state1)
        σ2_expected = MechMat.calculate_stress(me, ϵ2) + 2 * chain.G * (dev(ϵ2) - state1.ϵv[1])
        @test σ2_frozen ≈ σ2_expected
        σ2_true, _, _ = material_response(m, ϵ2, state1, 0.5)
        @test !(σ2_true ≈ σ2_frozen)

        # Reduced stress state: previously unsupported, now works automatically
        # via the generic fallback (GeneralizedMaxwell has no dedicated reduced
        # method, only the full-dimensional one used above).
        rss = ReducedStressState(PlaneStress(), m)
        ϵ1_red = SymmetricTensor{2,2}((0.01, 0.0, 0.0))
        state0_red = initial_material_state(rss)
        σ1_red, _, state1_red, _ = material_response(rss, ϵ1_red, state0_red, 0.5)
        @test stress_from_state(rss, ϵ1_red, state1_red) ≈ σ1_red
    end

    @testset "RotatedMaterial" begin
        e = LinearElastic(E=210.e3, ν=0.3)
        m_plastic = Plastic(elastic=e, yield=100.0, isotropic=Voce(Hiso=10.e3, κ∞=200.0), kinematic=ArmstrongFrederick(Hkin=1.e4, β∞=150.0))
        r = 2 * π * rand(Vec{3})
        rm = RotatedMaterial(m_plastic, r)

        state0 = initial_material_state(rm)
        ϵ_global1 = SymmetricTensor{2,3}((i, j) -> (i, j) == (1, 1) ? 0.01 : 0.0)
        _, _, state1 = material_response(rm, ϵ_global1, state0, nothing)

        ϵ_global2 = ϵ_global1 + SymmetricTensor{2,3}((i, j) -> (i, j) == (1, 1) ? 0.001 : 0.0)
        σ_current = stress_from_state(rm, ϵ_global2, state1)

        θ = norm(r)
        ϵ_local2 = rotate(ϵ_global2, r, -θ)
        σ_local_expected = MechMat.calculate_stress(e, ϵ_local2 - state1.ϵp)
        σ_expected = rotate(σ_local_expected, r, θ)
        @test σ_current ≈ σ_expected

        # Stateless wrapped material: verify no dispatch ambiguity and correct rotation
        m_el = LinearElastic{:cubicsymmetry}(C1111=1 + rand(), C1122=1 + rand(), C1212=1 + rand())
        rm_el = RotatedMaterial(m_el, r)
        ϵ = rand(SymmetricTensor{2,3})
        σ_rm_el = stress_from_state(rm_el, ϵ, initial_material_state(rm_el))
        σ_local_el = stress_from_state(m_el, rotate(ϵ, r, -θ), initial_material_state(m_el))
        @test σ_rm_el ≈ rotate(σ_local_el, r, θ)

        # Reduced stress state wrapping a rotated, stateful material: previously
        # unsupported, now works automatically via the generic fallback.
        rss = ReducedStressState(PlaneStress(), rm)
        ϵ1_red = SymmetricTensor{2,2}((0.01, 0.0, 0.0))
        state0_red = initial_material_state(rss)
        σ1_red, _, state1_red, _ = material_response(rss, ϵ1_red, state0_red, nothing)
        @test stress_from_state(rss, ϵ1_red, state1_red) ≈ σ1_red
    end

    @testset "HyperElastic" begin
        models = (NeoHooke(G=1 + rand()), CompressibleNeoHooke(G=1 + rand(), K=10 + rand()), SaintVenant(LinearElastic(E=210.e3, ν=0.3)))
        for m in models
            state = initial_material_state(m)
            F = one(Tensor{2,3}) + rand(Tensor{2,3}) / 20
            P, _, _ = material_response(m, F, state)
            @test stress_from_state(m, F, state) ≈ P

            # Reduced stress state: this already works via the generic
            # `NoMaterialState` fallback, since MaterialModelsBase's stress-state
            # iteration machinery already supports finite-strain (Tensor{2,3})
            # reduced states generically.
            rss = ReducedStressState(PlaneStress(), m)
            F_red = one(Tensor{2,2}) + rand(Tensor{2,2}) / 20
            state_red = initial_material_state(rss)
            P_red, _, _, _ = material_response(rss, F_red, state_red)
            @test stress_from_state(rss, F_red, state_red) ≈ P_red
        end
    end

    @testset "FiniteStrainPlastic" begin
        E, ν, Y0 = 210.e3, 0.3, 100.0
        nh = CompressibleNeoHooke(G=convert_hooke_param(:G; E, ν), K=convert_hooke_param(:K; E, ν))
        m = FiniteStrainPlastic(elastic=nh, yield=Y0, isotropic=Voce(Hiso=10.e3, κ∞=200.0), kinematic=ArmstrongFrederick(Hkin=1.e4, β∞=150.0))

        # Load to a converged, plastically loaded state
        state0 = initial_material_state(m)
        F1 = Tensor{2,3}((i, j) -> i == j ? (i == 1 ? 1.02 : 1.0) : 0.0)
        P1, _, state1 = material_response(m, F1, state0, nothing)
        @test stress_from_state(m, F1, state1) ≈ P1
        @test state1.Fp != state0.Fp # sanity: this test only matters if plastic loading occurred

        # Frozen-state postprocessing: a different F should give the frozen-Fp
        # elastic response, NOT a fresh plastic correction.
        F2 = Tensor{2,3}((i, j) -> i == j ? (i == 1 ? 1.03 : 1.0) : 0.0)
        σ2_frozen = stress_from_state(m, F2, state1)
        σ2_true, _, state2_true = material_response(m, F2, state1, nothing)
        @test !(σ2_true ≈ σ2_frozen) # material_response would further evolve plastically
        @test state2_true.Fp != state1.Fp

        # Reduced stress state, via the FrozenStressMaterial + MMB stress-state iteration
        rss = ReducedStressState(PlaneStress(), m)
        F1_red = Tensor{2,2}((1.02, 0.0, 0.0, 1.0))
        state0_red = initial_material_state(rss)
        P1_red, _, state1_red, _ = material_response(rss, F1_red, state0_red, nothing)
        @test stress_from_state(rss, F1_red, state1_red) ≈ P1_red
    end
end
