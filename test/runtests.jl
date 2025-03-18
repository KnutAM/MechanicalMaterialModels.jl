using Tensors, ForwardDiff, FiniteDiff, LinearAlgebra

using MechanicalMaterialModels
import MechanicalMaterialModels as MechMat

using MaterialModelsBase
import MaterialModelsBase as MMB

using Test

# Include files used for testing
include("utilities4testing.jl")

# Test of utility features
include("test_crystals.jl")

# Test material behaviors
include("test_elastic.jl")
include("test_hyperelastic.jl")
include("test_yield.jl")
include("test_hardening.jl")
include("test_overstress.jl")
include("test_plastic.jl")
include("test_viscoplastic.jl")
include("test_differentiate.jl")

# Test wrappers
include("test_rotated.jl")