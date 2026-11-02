@testset "FiniteStrainPlastic" begin
    # Test similarity to small strain implementation, using very high elastic stiffness / yield limit should provide similar results.
    E = 1000e3 # MPa
    ν = 0.3    # -
    Y0 = 10.0  # MPa
    Hkin = 300e3 * (1 + rand())
    β∞ = 20 * (1 + rand())
    κ∞ = 20 * (1 + rand())
    Hiso = 300e3 * (1 + rand())

    nh = CompressibleNeoHooke(;G = convert_hooke_param(:G; E, ν), K = convert_hooke_param(:K; E, ν))
    el = LinearElastic(;E, ν)
    kinematic = ArmstrongFrederick(;Hkin, β∞)
    isotropic = Voce(;Hiso, κ∞)

    ϵmax = 0.0001
    nstp = 100
    ϵ11_vec = append!(collect(range(0, ϵmax, nstp)), collect(range(ϵmax, -ϵmax, 2 * nstp)[2:end]))
    F11_vec = ϵ11_vec .+ 1

    m_fs = FiniteStrainPlastic(;elastic=nh, yield=Y0, kinematic=(kinematic,), isotropic=(isotropic,))
    m_ss = Plastic(;elastic=el, yield=Y0, kinematic=(kinematic,), isotropic=(isotropic,))

    σ11, σ, dσdϵ, state_ss, ϵ_full = run_normal(m_ss, ϵ11_vec; Δt = NaN, stress_state=UniaxialNormalStress(), ϵ_full = zero(SymmetricTensor{2,3}))
    P11, P, dPdF, state_fs, F_full = run_normal(m_fs, F11_vec; Δt = NaN, stress_state=UniaxialNormalStress(), ϵ_full = one(Tensor{2,3}))
    @test norm(σ11 - P11)/norm(P11) < 1e-4
end

@testset "FiniteStrainPlastic anisotropic non-symmetric Mandel stress" begin
    # `SaintVenant` wrapping an anisotropic `LinearElastic` produces a Mandel stress
    # M = Ce⋅Se that is genuinely non-symmetric (Se = C:E is not coaxial with Ce, unlike
    # for an isotropic elastic law). The VonMises flow direction used for Fp/backstress
    # evolution must then be computed consistently with the true (non-symmetric-aware)
    # gradient of the scalar effective stress.
    el_aniso = LinearElastic{:cubicsymmetry}(C1111=1200.0, C1122=400.0, C1212=300.0)
    elastic = SaintVenant(el_aniso)

    F = Tensor{2,3}((1.05, 0.08, -0.03, 0.02, 0.97, 0.06, -0.04, 0.01, 1.03))
    M = MechMat.mandel_stress(elastic, F)
    @test !(M ≈ M')  # confirms this case actually exercises the non-symmetric path

    yc = VonMises(100.0)
    ν = MechMat.effective_stress_gradient(yc, M)
    ν_ad = gradient(s -> MechMat.effective_stress(yc, s), M)
    @test ν ≈ ν_ad
end
