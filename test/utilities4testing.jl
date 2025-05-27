function run_shear(m, γmax, numsteps, Δt = nothing; TT=SymmetricTensor)
    state = initial_material_state(m)
    Δϵ = TT{2,3}((i,j)-> i==2 && j==1 ? γmax/numsteps : zero(γmax))
    s21 = zeros(numsteps+1)
    local σ, dσdϵ, ϵ
    for i = 1:numsteps
        ϵ = i*Δϵ
        σ, dσdϵ, state = material_response(m, ϵ, state, Δt)
        s21[i+1] = σ[2,1]
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

function obtain_numerical_material_derivative!(deriv, m, ϵ, old, Δt)
    p = tovector(m)
    # σ, dσdϵ, state = material_response(m, ϵ, old, Δt)
    dmdz = (σ = (ϵ = deriv.dσdϵ, s = deriv.dσdⁿs, p = deriv.dσdp),
            s = (ϵ = deriv.dsdϵ, s = deriv.dsdⁿs, p = deriv.dsdp))
    _tovector(t::AbstractTensor) = tomandel(t)
    _tovector(s::AbstractMaterialState) = tovector(s)
    funs = Tuple((ϵ=evec -> _tovector(material_response(m, frommandel(typeof(ϵ),evec), old, Δt)[i]), # f(ϵ)
                  s=svec -> _tovector(material_response(m, ϵ, MMB.fromvector(svec, old), Δt)[i]),    # f(s)
                  p=pvec -> _tovector(material_response(fromvector(pvec, m), ϵ, old, Δt)[i]))    # f(p)
                for i in [1,3])
    x0 = (ϵ = tomandel(ϵ), s = tovector(old), p = p)
    for (i, m) in enumerate((:σ,:s))
        tmp1 = getfield(dmdz,m)
        for z in fieldnames(typeof(x0))
            if length(getfield(x0, z)) > 0
                tmp=getfield(tmp1,z)
                fun(x) = getfield(funs[i], z)(x)
                tmp .= FiniteDiff.finite_difference_jacobian(fun, getfield(x0, z))
            end
        end
    end
end

function obtain_numerical_material_derivative!(deriv, m::LinearElastic, ϵ, old, Δt)
    # Set all derivatives to zero
    for field in fieldnames(typeof(deriv))
        fill!(getfield(deriv, field), 0)
    end
    dmdz = (σ = (ϵ = deriv.dσdϵ, s = deriv.dσdⁿs, p = deriv.dσdp),
            s = (ϵ = deriv.dsdϵ, s = deriv.dsdⁿs, p = deriv.dsdp))
    # σ, dσdϵ, state = material_response(m, ϵ, old, Δt)
    funs = (ϵ=evec -> tomandel(material_response(m, frommandel(typeof(ϵ),evec), old, Δt)[1]),          # f(ϵ)
            p=pvec -> tomandel(material_response(fromvector(pvec,m), ϵ, old, Δt)[1])) # f(p)
    x0 = (ϵ=tomandel(ϵ), p=tovector(m))
    
    dσdz = getfield(dmdz,:σ)
    for z in fieldnames(typeof(x0))
        tmp=getfield(dσdz,z)
        fun(x) = getfield(funs, z)(x)
        tmp .= FiniteDiff.finite_difference_jacobian(fun, getfield(x0, z))
    end
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
