@testset "OverstressFunction" begin
    @testset "RateIndependent" begin
        of = RateIndependent()
        @test get_vector_length(of) == 0
        @test of == fromvector(rand(10), of)
        
        @test contains(show_as_string(of), "Rate independent response")
    end

    @testset "NortonOverstress" begin
        tstar, nexp = rand(2)
        of = NortonOverstress(;tstar, nexp)
        @test get_vector_length(of) == 2
        
        MatTest.test_vectorconversion(Float64, of)
        MatTest.test_vectorconversion(Float32, of)
        MatTest.test_vectorconversion(MatTest.DualT{Float32}, of)

        Φ, σy = rand(2)
        # Scaling of Φ and σy doesn't change values
        @test MechMat.overstress_function(of, Φ, σy) ≈ MechMat.overstress_function(of, π*Φ, π*σy)
        # Exponent correctly implemented
        @test MechMat.overstress_function(of, π*Φ, σy) ≈ (π^of.nexp)*MechMat.overstress_function(of, Φ, σy)

        @test contains(show_as_string(of), "Norton overstress with")
    end
end