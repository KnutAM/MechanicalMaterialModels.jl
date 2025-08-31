function run_shear(m, γmax::Number, numsteps::Int, Δt = nothing; TT=SymmetricTensor)
    γv = range(0, γmax, numsteps + 1)
    return run_shear(m, γv, Δt; TT)
end

function run_shear(m, γv::AbstractVector, Δt = nothing; TT = SymmetricTensor)
    state = initial_material_state(m)
    s21 = zeros(length(γv))
    local σ, dσdϵ, ϵ
    for (k, ϵ21) in enumerate(γv[2:end])
        ϵ = TT{2,3}((i,j)-> i==2 && j==1 ? ϵ21 : zero(ϵ21))
        σ, dσdϵ, state = material_response(m, ϵ, state, Δt)
        s21[k+1] = σ[2,1]
    end
    return s21, σ, dσdϵ, state, ϵ
end

run_normal(m, ϵmax::Number, numsteps; kwargs...) = run_normal(m, Vector(range(0,ϵmax;length=numsteps+1)); kwargs...)

function run_normal(m, ϵ11_vec; Δt = NaN, stress_state=UniaxialNormalStress(), ϵ_full = zero(SymmetricTensor{2,3}))
    s11 = zeros(length(ϵ11_vec))
    #local ϵ_full
    local σ, dσdϵ, state
    old_state = initial_material_state(m)
    cache = allocate_material_cache(m)
    for (i, ϵ11) = enumerate(ϵ11_vec)
        ϵ = typeof(ϵ_full)((i,j)-> i==j==1 ? ϵ11 : ϵ_full[i,j])
        σ, dσdϵ, state, ϵ_full = material_response(stress_state, m, ϵ, old_state, Δt, cache)
        old_state = state
        s11[i] = σ[1,1]
    end
    return s11, σ, dσdϵ, state, ϵ_full
end

function test_conversion(mbase)
    v0 = [rand() for _ in 1:MMB.get_num_params(mbase)]
    v1 = similar(v0)
    m = fromvector(v0, mbase)
    tovector!(v1, m)
    @test v0 ≈ v1
end

function show_as_string(value, mime=MIME"text/plain"())
    io = IOBuffer()
    show(IOContext(io), mime, value)
    return String(take!(io))
end
