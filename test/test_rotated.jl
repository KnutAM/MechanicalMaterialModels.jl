@testset "RotatedMaterial" begin
    m_el = LinearElastic{:cubicsymmetry}(;C1111 = 1 + rand(), C1122 = 1 + rand(), C1212 = 1 + rand())
    r = 2 * π * rand(Vec{3})

    m = RotatedMaterial(m_el, r)
    D_rot = rotate(m_el.C, r, norm(r))

    ϵ = rand(SymmetricTensor{2,3})

    σ, dσdϵ, state = material_response(m, ϵ, initial_material_state(m))

    @test σ ≈ D_rot ⊡ ϵ
    @test dσdϵ ≈ D_rot

    @testset "zero rotation vector" begin
        m0 = RotatedMaterial(m_el, zero(Vec{3}))
        σ0, dσ0dϵ, _ = material_response(m0, ϵ, initial_material_state(m0))
        σ_direct, dσdϵ_direct, _ = material_response(m_el, ϵ, initial_material_state(m_el))

        @test !any(isnan, σ0)
        @test !any(isnan, dσ0dϵ)
        @test σ0 ≈ σ_direct
        @test dσ0dϵ ≈ dσdϵ_direct
    end

    @testset "tiny (near-zero) rotation vector" begin
        r_tiny = Vec{3}((1.0e-9, 0.0, 0.0))
        m_tiny = RotatedMaterial(m_el, r_tiny)
        σ_tiny, dσ_tinydϵ, _ = material_response(m_tiny, ϵ, initial_material_state(m_tiny))
        σ_direct, _, _ = material_response(m_el, ϵ, initial_material_state(m_el))

        @test !any(isnan, σ_tiny)
        @test !any(isnan, dσ_tinydϵ)
        # A vanishingly small rotation should give (almost) the unrotated response
        @test σ_tiny ≈ σ_direct atol = 1.0e-5

        # Extreme underflow: norm(r) itself underflows to exactly 0.0 in Float64,
        # so this exercises the Taylor-series branch, not just the exact-zero path.
        r_underflow = Vec{3}((1.0e-200, 0.0, 0.0))
        @test iszero(norm(r_underflow))
        m_underflow = RotatedMaterial(m_el, r_underflow)
        σ_uf, dσ_ufdϵ, _ = material_response(m_underflow, ϵ, initial_material_state(m_underflow))
        @test !any(isnan, σ_uf)
        @test !any(isnan, dσ_ufdϵ)
        @test σ_uf ≈ σ_direct atol = 1.0e-5
        @test stress_from_state(m_underflow, ϵ, initial_material_state(m_underflow)) ≈ σ_uf
    end

    @testset "non-Float64 and integer rotation vectors" begin
        # Regression: a naive `oftype(θsq, threshold)` guard throws for integer
        # rotation vectors and silently rounds the threshold to zero for Float16.
        for r_special in (zero(Vec{3, Int}), zero(Vec{3, Float16}), Vec{3}((0, 0, 1)))
            m_special = RotatedMaterial(m_el, r_special)
            σ_special, _, _ = material_response(m_special, ϵ, initial_material_state(m_special))
            @test !any(isnan, σ_special)
        end
    end

    @testset "orientation derivative through zero rotation" begin
        # The rotation tensor must remain differentiable with respect to the
        # rotation vector itself, including exactly at r = 0 (a naive
        # `iszero(r) && return unrotated response` guard would give the wrong,
        # zero derivative here).
        m_aniso = LinearElastic{:cubicsymmetry}(C1111=4.0, C1122=1.0, C1212=0.7)
        ϵ_aniso = SymmetricTensor{2,3}((1.0, 0.2, 0.0, 0.3, 0.0, 0.0))
        f(t) = material_response(RotatedMaterial(m_aniso, Vec{3}((zero(t), zero(t), t))),
            ϵ_aniso, initial_material_state(m_aniso))[1][1, 2]

        d_ad = ForwardDiff.derivative(f, 0.0)
        d_fd = (f(1.0e-6) - f(-1.0e-6)) / 2.0e-6

        @test !isnan(d_ad)
        @test d_ad ≈ d_fd rtol = 1.0e-4
        @test abs(d_ad) > 1.0e-3 # sanity: this test only matters if the derivative is nonzero
    end

    @testset "finite-strain (non-symmetric) tensors" begin
        F = one(Tensor{2,3}) + rand(Tensor{2,3}) / 20
        rm0 = RotatedMaterial(CompressibleNeoHooke(G=80.0e3, K=160.0e3), zero(Vec{3}))
        P0, dPdF0, _ = material_response(rm0, F, initial_material_state(rm0))
        P_direct, dPdF_direct, _ = material_response(rm0.material, F, initial_material_state(rm0.material))

        @test !any(isnan, P0)
        @test !any(isnan, dPdF0)
        @test P0 ≈ P_direct
        @test dPdF0 ≈ dPdF_direct
    end
end
