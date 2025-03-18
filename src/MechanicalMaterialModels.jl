module MechanicalMaterialModels
using Tensors, StaticArrays, ForwardDiff, LinearAlgebra

import Base: @kwdef  # Note exported for Julia version < 1.9

using MaterialModelsBase, Newton
import MaterialModelsBase as MMB

# Utils etc.
include("utils.jl")

# Wrappers
include("RotatedMaterial.jl")
export RotatedMaterial

# Hardening laws etc. 
include("plasticity_components/YieldCriteria.jl")
export VonMises, DruckerPrager

include("plasticity_components/IsotropicHardening.jl")
export Voce, Swift

include("plasticity_components/KinematicHardening.jl")
export ArmstrongFrederick, Delobelle, OhnoWang

include("plasticity_components/OverstressFunctions.jl")
export RateIndependent, NortonOverstress

include("plasticity_components/Crystallography.jl")
export BCC, BCC12, FCC, GenericCrystallography

# Different material classes
include("Elastic.jl")
export LinearElastic

include("Plastic.jl")
include("ExtraOutputs.jl")
include("PlasticDifferentiate.jl")
export Plastic

include("hyper_elasticity/HyperElastic.jl")
include("hyper_elasticity/NeoHooke.jl")
export NeoHooke, CompressibleNeoHooke

include("FiniteStrainPlastic.jl")

end
