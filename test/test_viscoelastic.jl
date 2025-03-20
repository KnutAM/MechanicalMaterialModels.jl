@testset "Viscoelasticity" begin
    G0 = 10 * (1 + rand())
    K0 = 20 * (1 + rand())
    Gs = ntuple(_ -> G0 * (0.5 + rand()), 5)         # Random chain stiffness around G0
    ts = ntuple(i -> (1 + rand()) * 10.0 ^ (i - 2), 5) # Relaxation times ~ 0.01, 0.1, 1, 10, 100
    E0 = 9 * G0 * K0 / (3 * K0 + G0)
    ν0 = E0 / (2 * G0) - 1
    me = LinearElastic(;E = E0, ν = ν0)
    chains = map((G, t) -> Maxwell(;G, t), Gs, ts)

    m1 = GeneralizedMaxwell(me, chains[1])
    m2 = GeneralizedMaxwell(me, chains[1:2])
    m5 = GeneralizedMaxwell(me, chains)

    function instant_stiffness(m::GeneralizedMaxwell)
        I2 = one(SymmetricTensor{2,3})
        I4dev = one(SymmetricTensor{4,3}) - I2 ⊗ I2 / 3
        return m.base.C + sum(c -> 2 * c.G, m.chains) * I4dev
    end
    relaxation_time(m::Maxwell) = m.η / m.G
    relaxation_time(m::GeneralizedMaxwell) = maximum(relaxation_time, m.chains)
    
    # Instantanious and long-term stiffness responses
    D0 = me.C
    @testset "GM (num = $(length(m.chains)))" for m in (m1, m2, m5)
        Δϵ = rand(SymmetricTensor{2,3}) / 100
        state = initial_material_state(m)
        # Instant stiffness
        σ, dσdϵ, _ = material_response(m, Δϵ, state, 1e-20)
        @test dσdϵ ≈ instant_stiffness(m)
        @test σ ≈ instant_stiffness(m) ⊡ Δϵ
        # Long-term secant stiffness
        σ, dσdϵ, _ = material_response(m, Δϵ, state, 1e20)
        @test dσdϵ ≈ D0
        @test σ ≈ D0 ⊡ Δϵ

        # Consider load stepping and check that we get the correct saturated stress level 
        T_max = 20 * relaxation_time(m)
        n_steps = 100
        t_old = 0
        state = initial_material_state(m)
        ϵ_end = rand(SymmetricTensor{2,3}) / 100
        for i in 1:n_steps
            ϵ = clamp(i/10, 0, 1) * ϵ_end
            t = T_max * (i / n_steps) ^ 2
            Δt = t - t_old
            σ, _, state = material_response(m, ϵ, state, Δt)
            t_old = t
        end
        @test σ ≈ D0 ⊡ ϵ_end
    end

    # Test some error paths and constructors 
    η = rand()
    c1 = Maxwell(;G = G0, η)
    c2 = Maxwell(;G = G0, t = η / G0)
    @test_throws ArgumentError Maxwell(;G = G0)
    @test_throws ArgumentError Maxwell(;G = G0, η, t = rand())

    @test c1.η == η
    @test c2.η ≈ η # Correct conversion
end
