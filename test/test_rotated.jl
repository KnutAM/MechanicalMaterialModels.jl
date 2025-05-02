@testset "RotatedMaterial" begin
    m_el = LinearElastic{:cubicsymmetry}(;C1111 = 1 + rand(), C1122 = 1 + rand(), C1212 = 1 + rand())
    r = 2 * π * rand(Vec{3})

    m = RotatedMaterial(m_el, r)
    D_rot = rotate(m_el.C, r, norm(r))
    
    ϵ = rand(SymmetricTensor{2,3})

    σ, dσdϵ, state = material_response(m, ϵ, initial_material_state(m))

    @test σ ≈ D_rot ⊡ ϵ
    @test dσdϵ ≈ D_rot
end
