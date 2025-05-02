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

    s11, σ, dσdϵ, state, ϵ = run_normal(m, Δϵ21, N; Δt = 1e-3)
end
