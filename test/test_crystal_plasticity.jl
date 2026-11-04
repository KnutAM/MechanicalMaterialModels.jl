@testset "CrystalPlasticity" begin
    m = CrystalPlasticity(;
        crystal = FCC(),
        elastic = LinearElastic(;E = 100.e3, ν = 0.3),
        yield = 50.0,
        q = 0.5, h0 = 10e3, h∞ = 1e3, ζ = 1.0,
        Hkin = 0.0, β∞ = 25.0,
        overstress = NortonOverstress(;tstar = 1.0, nexp = 2.0))
    s = initial_material_state(m)

    Δϵ21 = 0.01
    N = 100
    s21, σ, dσdϵ, state, ϵ = run_shear(m, Δϵ21, N, 1e-3)
    @test state.γ_acc >= 0
    @test issorted(s21)

    s11, σ, dσdϵ, state, ϵ = run_normal(m, Δϵ21, N; Δt = 1e-3)
    @test state.γ_acc >= 0

    @testset "overstress must be rate-dependent (bug 7)" begin
        # RateIndependent forces Φᵢ=0 on every slip system simultaneously once any
        # system yields, including inactive ones - this can drive an inactive
        # system to a spurious negative slip increment. Disallowed at construction
        # instead, requiring a viscoplastic (e.g. Norton) regularization, whose
        # residual is non-negative by construction for an inactive system.
        @test_throws MethodError CrystalPlasticity(;
            crystal = FCC(), elastic = LinearElastic(;E = 100.e3, ν = 0.3),
            yield = 50.0, q = 0.5, h0 = 10e3, h∞ = 1e3, ζ = 1.0,
            Hkin = 0.0, β∞ = 25.0, overstress = RateIndependent())
    end

    @testset "frozen slip direction vs. converged stress sign (bug 8)" begin
        # Two coupled slip systems (Schmid tensor coupling 0.8): the trial
        # (elastic-predictor) resolved shear stress on both systems is positive,
        # so both get ν_trial=+1 frozen for the whole Newton solve. Cross-system
        # elastic coupling flips the *converged* sign of the resolved shear on
        # system 2 during the correction, but the old code didn't check this and
        # silently converged to a positive Δγ for system 2 using the wrong
        # (stale) frozen direction, yielding negative plastic dissipation - a
        # thermodynamic violation. Should now be rejected instead of silently
        # accepted.
        θ = acos(0.8) / 2
        planes = MechMat.SVector(Vec{3}((0.0, 1.0, 0.0)), Vec{3}((-sin(θ), cos(θ), 0.0)))
        directions = MechMat.SVector(Vec{3}((1.0, 0.0, 0.0)), Vec{3}((cos(θ), sin(θ), 0.0)))
        crystal = GenericCrystallography(planes, directions)
        mc = CrystalPlasticity(; crystal, elastic = LinearElastic(E = 2.0, ν = 0.0), yield = 0.5,
            q = 0.0, h0 = 0.0, h∞ = 0.0, ζ = 1.0, Hkin = 0.0, β∞ = 1.0,
            overstress = NortonOverstress(tstar = 2.0, nexp = 1.0))
        ϵ = SymmetricTensor{2,3}([35/6 5.0 0.0; 5.0 -35/6 0.0; 0.0 0.0 0.0])
        old = initial_material_state(mc)
        @test_throws MMB.NoLocalConvergence material_response(mc, ϵ, old, 1.0)

        # Edge case: the frozen trial direction can itself be exactly zero
        # (ν_trial[i] = sign(0) = 0) when the trial resolved shear on a system
        # exactly cancels its backstress. A naive `ν_trial[i]*τ_conv[i] >= -tol`
        # check would then pass trivially regardless of the converged stress,
        # silently accepting an inconsistent (zero) flow direction for a system
        # that is, in fact, significantly active.
        σ_trial = MechMat.calculate_stress(mc.elastic, ϵ - old.ϵp)
        slip_dirs = MechMat.get_slip_directions(crystal)
        slip_planes = MechMat.get_slip_planes(crystal)
        τ_trial_red = map((β, s, p) -> s ⋅ σ_trial ⋅ p - β, old.β, slip_dirs, slip_planes)
        β2 = τ_trial_red[2] + old.β[2] # makes the adjusted τ_trial_red[2] exactly 0
        old_zero_trial = MechMat.CrystalPlasticityState(old.γ_acc, old.ϵp, old.κ, MechMat.SVector(old.β[1], β2))
        @test_throws MMB.NoLocalConvergence material_response(mc, ϵ, old_zero_trial, 1.0)

        # Edge case: the check must be scaled by the stress-like yield/backstress
        # quantities, not compared against `mat.tolerance` directly (which bounds
        # a dimensionless slip residual). Otherwise rescaling every stress-valued
        # parameter (E, yield, β∞) down lets the same inconsistent solution pass
        # silently, since the raw stress mismatch then falls below the
        # (unscaled) tolerance.
        scale = 1e-10
        mc_scaled = CrystalPlasticity(; crystal, elastic = LinearElastic(E = 2.0*scale, ν = 0.0), yield = 0.5*scale,
            q = 0.0, h0 = 0.0, h∞ = 0.0, ζ = 1.0, Hkin = 0.0, β∞ = 1.0*scale,
            overstress = NortonOverstress(tstar = 2.0, nexp = 1.0))
        old_scaled = initial_material_state(mc_scaled)
        @test_throws MMB.NoLocalConvergence material_response(mc_scaled, ϵ, old_scaled, 1.0)
    end
end
