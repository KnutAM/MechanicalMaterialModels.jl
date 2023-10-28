@testset "ViscoPlasticCombined" begin
    elast = LinearElastic(; E=210.e3, ν=0.3)
    σ_y0 = 200.
    kin = ArmstrongFrederick(;Hkin=10.e3, β∞=300.0)
    iso = Voce(;Hiso=5.e3, κ∞=400.0)
    
    m_rate_indep = Plastic(elastic=elast, yield=σ_y0, isotropic=iso, kinematic=kin)
    m_norton = Plastic(elastic=elast, yield=σ_y0, isotropic=iso, kinematic=kin, overstress=NortonOverstress(;nexp=2, tstar=1e-3))

    function run_normal_named(args...; kwargs...)
        s11, σ, dσdϵ, state, ϵ = run_normal(args...; kwargs...)
        return (s11=s11, σ=σ, dσdϵ=dσdϵ, state=state, ϵ=ϵ)
    end

    ϵ11 = 0.02
    numsteps = 20
    Δt=0.05

    res_rate_indep = run_normal_named(m_rate_indep, ϵ11, numsteps; Δt=Δt)

    # Rate dependent, stress should be higher
    res_norton = run_normal_named(m_norton, ϵ11, numsteps; Δt=Δt)
    
    @test !(isapprox(res_norton.s11, res_rate_indep.s11; rtol=1e-5))    # Using same tol as below to ensure sufficient diff. 
    @test all(res_norton.s11 .>= res_rate_indep.s11)

    res_norton_slow = run_normal_named(m_norton, ϵ11, numsteps; Δt=Δt*1e6)
    @test isapprox(res_norton_slow.s11, res_rate_indep.s11; rtol=1e-5)

    # Conversions
    function check_conversions(m)
        v = rand(get_num_params(m))
        mr = vector2material(v, m)
        vc = material2vector(mr)
        @test vc ≈ v 
        fill!(vc, rand())
        material2vector!(vc, mr)
        @test vc ≈ v
    end

    check_conversions(m_rate_indep)
    check_conversions(m_norton)
end