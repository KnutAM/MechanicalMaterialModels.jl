"""
    RotatedMaterial(m::AbstractMaterial, r::Vec) <: AbstractMaterial

A rotated material is wrapper around `m`, such that the material response is evaluated in a 
local coordinate system defined by the Rodrigues rotation vector, `r`. Specifically, the input 
strain is first rotated by `θ = |r|` around `r` before calling `material_response` with `m`. 
The resulting stress and stiffness are rotated back, `θ = -|r|` around `r`, before they are returned. 

!!! note State variables are in the local coordinate system
    The state variables are not modified, and are hence defined in the local coordinate system. 
    Care must therefore be taken to rotate the strain to the local coordinates when evaluating the 
    responses, and the output of those evaluations should be rotated back to the global coordinates.

"""
struct RotatedMaterial{M <: AbstractMaterial, RT <: Vec} <: AbstractMaterial
    material::M
    rotation::RT
end

for op in ( :initial_material_state, :allocate_material_cache,
            :get_tensorbase, :get_num_statevars,
            :allocate_differentiation_output)
    @eval @inline MMB.$op(rm::RotatedMaterial) = MMB.$op(rm.material)
end

function MMB.material_response(rm::RotatedMaterial, strain::AbstractTensor, args::Vararg{Any, N}) where {N}
    θ = norm(rm.rotation)
    strain_rot = rotate(strain, rm.rotation, θ)
    stress_rot, stiff_rot, state = MMB.material_response(m::Plastic, strain_rot::SymmetricTensor{2,3}, args...)
    stress = rotate(stress_rot, rm.rotation, -θ)
    stiff = rotate(stiff_rot, rm.rotation, -θ)
    return stress, stiff, state
end
