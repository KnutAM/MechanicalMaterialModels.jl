@testset "HyperElastic" begin
    models = (NeoHooke(;G=rand()), CompressibleNeoHooke(;G=rand(), K=rand()))
    F = rand(Tensor{2,3})
    for model in models
        state = initial_material_state(model)
        ∂P∂F_ad, P_ad = Tensors.hessian(x->MechMat.compute_potential(model, tdot(x)), F, :all)
        P, ∂P∂F, state = material_response(model, F, state)
        @test P ≈ P_ad
        @test ∂P∂F ≈ ∂P∂F_ad
        @test isa(state, NoMaterialState)
    end
end