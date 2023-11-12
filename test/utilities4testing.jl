function run_shear(m, γmax, numsteps)
    state = initial_material_state(m)
    Δϵ = SymmetricTensor{2,3}((i,j)-> i==2 && j==1 ? γmax/numsteps : zero(γmax))
    s21 = zeros(numsteps+1)
    local σ, dσdϵ, ϵ
    for i = 1:numsteps
        ϵ = i*Δϵ
        σ, dσdϵ, state = material_response(m, ϵ, state)
        s21[i+1] = σ[2,1]
    end
    return s21, σ, dσdϵ, state, ϵ
end

run_normal(m, ϵmax::Number, numsteps; kwargs...) = run_normal(m, Vector(range(0,ϵmax;length=numsteps+1)); kwargs...)

function run_normal(m, ϵ11_vec; Δt = NaN, stress_state=UniaxialNormalStress())
    s11 = zeros(length(ϵ11_vec))
    local ϵ_full = zero(SymmetricTensor{2,3})
    local σ, dσdϵ, state
    old_state = initial_material_state(m)
    cache = allocate_material_cache(m)
    for (i, ϵ11) = enumerate(ϵ11_vec)
        ϵ = SymmetricTensor{2,3}((i,j)-> i==j==1 ? ϵ11 : ϵ_full[i,j])
        σ, dσdϵ, state, ϵ_full = material_response(stress_state, m, ϵ, old_state, Δt, cache)
        old_state = state
        s11[i] = σ[1,1]
    end
    return s11, σ, dσdϵ, state, ϵ_full
end

function obtain_numerical_material_derivative!(deriv, m, ϵ, old, Δt)
    p = material2vector(m)
    # σ, dσdϵ, state = material_response(m, ϵ, old, Δt)
    dmdz = (σ = (ϵ = deriv.dσdϵ, s = deriv.dσdⁿs, p = deriv.dσdp),
            s = (ϵ = deriv.dsdϵ, s = deriv.dsdⁿs, p = deriv.dsdp))
    funs = Tuple((ϵ=evec -> tomandel(material_response(m, frommandel(typeof(ϵ),evec), old, Δt)[i]),          # f(ϵ)
                  s=svec -> tomandel(material_response(m, ϵ, frommandel(typeof(old),svec), Δt)[i]),          # f(s)
                  p=pvec -> tomandel(material_response(MMB.vector2material(pvec,m), ϵ, old, Δt)[i])) # f(p)
                for i in [1,3])
    x0 = (ϵ=tomandel(ϵ), s=tomandel!(zeros(get_num_statevars(m)),old), p=p)
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
            p=pvec -> tomandel(material_response(vector2material(pvec,m), ϵ, old, Δt)[1])) # f(p)
    x0 = (ϵ=tomandel(ϵ), p=material2vector(m))
    
    dσdz = getfield(dmdz,:σ)
    for z in fieldnames(typeof(x0))
        tmp=getfield(dσdz,z)
        fun(x) = getfield(funs, z)(x)
        tmp .= FiniteDiff.finite_difference_jacobian(fun, getfield(x0, z))
    end
end
