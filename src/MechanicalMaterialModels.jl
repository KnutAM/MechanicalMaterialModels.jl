module MechanicalMaterialModels
using Tensors, StaticArrays, ForwardDiff, LinearAlgebra

using MaterialModelsBase, Newton
import MaterialModelsBase as MMB

# Utils etc.
include("utils.jl")

# Hardening laws etc. 
include("ModelComponents/YieldCriteria.jl")
export VonMises, DruckerPrager

include("ModelComponents/IsotropicHardening.jl")
export Voce, Swift

include("ModelComponents/KinematicHardening.jl")
export ArmstrongFrederick, Delobelle, OhnoWang

include("ModelComponents/OverstressFunctions.jl")
export RateIndependent, NortonOverstress

# Different material classes
include("Elastic.jl")
export LinearElastic

include("Plastic.jl")
include("ExtraOutputs.jl")
include("PlasticDifferentiate.jl")
export Plastic

end
