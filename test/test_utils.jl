@testset "isotropic_hooke_conversions" begin
    E = 1 + rand()          # ∈ [1, 2]
    ν = 0.05 + 0.4 * rand() # ∈ [0.05, 0.45]
    K = convert_hooke_param(:K; E, ν)
    G = convert_hooke_param(:G; E, ν)
    @test K ≈ convert_hooke_param(:K; E, G)
    @test K ≈ convert_hooke_param(:K; G, ν)
    @test K ≈ convert_hooke_param(:K; G, E) # Reverse order
    @test K ≈ convert_hooke_param(:K; ν, E) # Reverse order

    @test G ≈ convert_hooke_param(:G; E, K)
    @test G ≈ convert_hooke_param(:G; ν, K)

    @test E ≈ convert_hooke_param(:E; ν, K)
    @test E ≈ convert_hooke_param(:E; ν, G)
    @test E ≈ convert_hooke_param(:E; G, K)
    
    @test ν ≈ convert_hooke_param(:ν; E, K)
    @test ν ≈ convert_hooke_param(:ν; E, G)
    @test ν ≈ convert_hooke_param(:ν; G, K)
     
end