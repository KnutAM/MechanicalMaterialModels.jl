# Examples

## Crystal plasticity
```@example crystalplasticity
using MaterialModelsBase, Tensors
using MechanicalMaterialModels
import CairoMakie as Plt

function simulate_uniaxial(m::AbstractMaterial, ϵ11_history, time_history)
    state = initial_material_state(m)
    cache = allocate_material_cache(m)
    stress_state = UniaxialStress()
    σ11_history = zero(ϵ11_history)
    for i in eachindex(ϵ11_history, time_history)[2:end]
        Δt = time_history[i] - time_history[i-1]
        ϵ = SymmetricTensor{2,1}((ϵ11_history[i],))
        σ, dσdϵ, state = material_response(stress_state, m, ϵ, state, Δt, cache)
        σ11_history[i] = σ[1,1]
    end
    return σ11_history, state
end

material = CrystalPlasticity(;
    crystal = FCC(),
    elastic = LinearElastic(;E = 100.e3, ν = 0.3),
    yield = 50.0,
    q = 0.5, h0 = 100e3, h∞ = 1e3, ζ = 1e3,
    Hkin = 0.0, β∞ = 25.0, 
    overstress = NortonOverstress(;tstar = 1.0, nexp = 2.0))

ϵ11_history  = collect(range(0, 0.01; length=100))  # Ramp to 1 %
time_history = collect(range(0, 100; length=length(ϵ11_history)))   # Constant time step
σ11_history, state = simulate_uniaxial(material, ϵ11_history, time_history)

fig = Plt.Figure()
ax = Plt.Axis(fig[1,1]; xlabel = "strain [%]", ylabel = "stress [MPa]")
Plt.lines!(ax, ϵ11_history * 100, σ11_history)
fig
```
