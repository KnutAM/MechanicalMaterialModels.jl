@testset "LinearElastic" begin
    # Test each constructor
    C4 = rand(SymmetricTensor{4,3})
    m0 = LinearElastic(C4)
    v1 = [1.0, 0.3]
    m1 = LinearElastic(E=v1[1]; ν=v1[2])
    v2 = [1.0, 0.9, 0.5]
    m2 = LinearElastic{:cubicsymmetry}(C1111=v2[1], C1122=v2[2], C1212=v2[3])

    # Test that construction via vector gives the same material
    @test m1 == fromvector(v1, m1)
    @test m2 == fromvector(v2, m2)

    # Test get stress function
    E = 210.e3; ν = 0.3
    G = E/(2*(1+ν)); K = E/(3*(1-2*ν))
    m = LinearElastic{:isotropic}(E=E, ν=ν)
    ϵ = rand(SymmetricTensor{2,3})
    σ_verify = 2G*dev(ϵ) + 3K*vol(ϵ)
    @test MechMat.calculate_stress(m, ϵ) ≈ σ_verify

    # Test material material_response
    σ, dσdϵ, state = material_response(m, ϵ, initial_material_state(m), 1.0, allocate_material_cache(m))
    @test σ ≈ σ_verify
    σfun(ϵ) = material_response(m, ϵ, initial_material_state(m))[1]
    @test Tensors.gradient(σfun, ϵ) ≈ dσdϵ

    # Show methods 
    @test contains(show_as_string(LinearElastic(E=2.0, ν=0.3)), "Isotropic")
    @test contains(show_as_string(LinearElastic{:general}(m.C)), "Fully anisotropic")
    m_cubic = LinearElastic{:cubicsymmetry}(; C1111=1., C1122=2., C1212=3.)
    @test contains(show_as_string(m_cubic), "Cubic symmetric")
    m_custom = LinearElastic{Float64,:custom,5}(m.C, rand(5))
    @test contains(show_as_string(m_custom), "custom")

    # Regression test: cubic-symmetry shear response must match an equivalent
    # isotropic material (shear stress should not be doubled)
    G_shear = 400.0
    m_cubic_shear = LinearElastic{:cubicsymmetry}(C1111=1200.0, C1122=400.0, C1212=G_shear)
    m_iso_shear = LinearElastic{:isotropic}(E=2*G_shear*(1+0.0), ν=0.0)
    ϵ12 = SymmetricTensor{2,3}((i,j) -> (i,j) in ((1,2),(2,1)) ? 0.01 : 0.0)
    σ_cubic = MechMat.calculate_stress(m_cubic_shear, ϵ12)
    σ_iso = MechMat.calculate_stress(m_iso_shear, ϵ12)
    @test σ_cubic[1,2] ≈ σ_iso[1,2]
    @test σ_cubic[1,2] ≈ 8.0

    # Regression test: the stiffness tensor must match the 9-component Voigt
    # matrix documented in the `LinearElastic{:cubicsymmetry}` docstring,
    # including the off-diagonal minor-symmetry couplings
    C1111, C1122, C1212 = v2
    m_doc = LinearElastic{:cubicsymmetry}(; C1111, C1122, C1212)
    voigt_doc = [
        C1111    C1122    C1122    0        0        0        0        0        0
        C1122    C1111    C1122    0        0        0        0        0        0
        C1122    C1122    C1111    0        0        0        0        0        0
        0        0        0        C1212    0        0        C1212    0        0
        0        0        0        0        C1212    0        0        C1212    0
        0        0        0        0        0        C1212    0        0        C1212
        0        0        0        C1212    0        0        C1212    0        0
        0        0        0        0        C1212    0        0        C1212    0
        0        0        0        0        0        C1212    0        0        C1212
    ]
    @test tovoigt(Tensor{4,3}(m_doc.C)) ≈ voigt_doc
end
