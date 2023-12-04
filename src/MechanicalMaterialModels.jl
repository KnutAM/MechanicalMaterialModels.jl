module MechanicalMaterialModels
using Tensors, StaticArrays, ForwardDiff, LinearAlgebra

import Base: @kwdef  # Note exported for Julia version < 1.9

using MaterialModelsBase, Newton
import MaterialModelsBase as MMB

# Utils etc.
include("utils.jl")

# Hardening laws etc. 
include("plasticity_components/YieldCriteria.jl")
export VonMises, DruckerPrager

include("plasticity_components/IsotropicHardening.jl")
export Voce, Swift

include("plasticity_components/KinematicHardening.jl")
export ArmstrongFrederick, Delobelle, OhnoWang

include("plasticity_components/OverstressFunctions.jl")
export RateIndependent, NortonOverstress

# Different material classes
include("Elastic.jl")
export LinearElastic

include("Plastic.jl")
include("ExtraOutputs.jl")
include("PlasticDifferentiate.jl")
export Plastic

end
