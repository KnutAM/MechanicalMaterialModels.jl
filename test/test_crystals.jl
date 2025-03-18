@testset "Crystallography.jl" begin
    # Could move this functionality into package to allow easy construction based on miller indices. 
    # However, should be noted that this code only works for cubic crystals, and not e.g. hexagonal. 
    function slip_planes(planes::NTuple{3, Int})
        all_planes = Set{NTuple{3, Int}}()
        remind(i,j) = (s = i+j; s == 3 ? 3 : (s == 4 ? 2 : 1))
        for i in 1:3
            for j in 1:3
                i == j && continue
                k = remind(i,j)
                p = (planes[i], planes[j], planes[k])
                for sign_ind in 1:3
                    v = ntuple(l -> p[l] * (l == sign_ind ? -1 : 1), 3)
                    (-).(v) ∉ all_planes && push!(all_planes, v)
                end
                (-).(p) ∉ all_planes && push!(all_planes, p)
            end
        end
        return sort(collect(all_planes))
    end
    function slip_directions(plane::NTuple{3, Int}, directions::NTuple{3, Int})
        all_dirs = Set{NTuple{3, Int}}()
        tupledot(a, b) = sum(a .* b)
        remind(i,j) = (s = i+j; s == 3 ? 3 : (s == 4 ? 2 : 1))
        for i in 1:3
            for j in 1:3
                i == j && continue
                k = remind(i, j)
                d = (directions[i], directions[j], directions[k])
                for sign_ind in 1:3
                    v = ntuple(l -> d[l] * (l == sign_ind ? -1 : 1), 3)
                    if tupledot(v, plane) == 0 # Perpendicular
                        (-).(v) ∉ all_dirs && push!(all_dirs, v)
                    end
                end
                if tupledot(d, plane) == 0 # Perpendicular
                    (-).(d) ∉ all_dirs && push!(all_dirs, d)
                end
            end
        end
        return sort(collect(all_dirs))
    end
    function get_slip_systems(;plane_set, direction_set)
        planes = slip_planes(plane_set)
        system_planes = similar(planes, 0)
        system_directions = similar(planes, 0)
        
        for plane in planes
            for direction in slip_directions(plane, direction_set)
                push!(system_planes, plane)
                push!(system_directions, direction)
            end
        end
        return map(normalize ∘ Vec, system_planes), map(normalize ∘ Vec, system_directions)
    end

    gen_c_2d = GenericCrystallography(2 * π * rand(3)...)
    @testset "Generic tests for $(typeof(c))" for c in (BCC(), BCC12(), FCC(), gen_c_2d)
        p = MechanicalMaterialModels.get_slip_planes(c)
        d = MechanicalMaterialModels.get_slip_directions(c)
        @test all(v -> norm(v) ≈ 1, p)
        @test all(v -> norm(v) ≈ 1, d)
        @test all(((x, y),) -> abs(x ⋅ y) < 1e-8, zip(p, d))
    end

    function check_same_systems(c; p_test, d_test)
        ps = MechanicalMaterialModels.get_slip_planes(c)
        ds = MechanicalMaterialModels.get_slip_directions(c)
        for (p, d) in zip(ps, ds)
            c = 0
            for (pt, dt) in zip(p_test, d_test)
                if (pt ≈ p || pt ≈ -p) && (dt ≈ d || dt ≈ -d)
                    c += 1
                end
            end
            @test c == 1
        end
    end

    @testset "FCC" begin
        p_test, d_test = get_slip_systems(;plane_set = (1,1,1), direction_set = (1,1,0))
        @test length(p_test) == length(d_test) == 12
        check_same_systems(FCC(); p_test, d_test)
    end

    @testset "BCC" begin
        p_test, d_test = get_slip_systems(;plane_set = (1,1,0), direction_set = (1,1,1))
        @test length(p_test) == length(d_test) == 12
        p2, d2 = get_slip_systems(;plane_set = (1,2,3), direction_set = (1,1,1))
        @test length(p2) == length(d2) == 24
        p3, d3 = get_slip_systems(;plane_set = (1,1,2), direction_set = (1,1,1))
        @test length(p3) == length(d3) == 12
        append!(p_test, p2)
        append!(p_test, p3)
        append!(d_test, d2)
        append!(d_test, d3)

        check_same_systems(BCC(); p_test, d_test)
    end

    @testset "BCC12" begin
        p_test, d_test = get_slip_systems(;plane_set = (1,1,0), direction_set = (1,1,1))
        @test length(p_test) == length(d_test) == 12
        check_same_systems(BCC12(); p_test, d_test)
    end
end
