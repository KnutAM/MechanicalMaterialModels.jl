using Tensors, ForwardDiff, FiniteDiff, LinearAlgebra

using MechanicalMaterialModels
import MechanicalMaterialModels as MechMat
import MaterialModelsTesting as MatTest

using MaterialModelsBase
import MaterialModelsBase as MMB

using Test

# Include files used for testing
include("utilities4testing.jl")

# Test of utility features
include("test_crystals.jl")
include("test_utils.jl")
include("test_hyperelastic.jl")

# Test modeling components
include("test_elastic.jl")
include("test_yield.jl")
include("test_hardening.jl")
include("test_overstress.jl")

# Test small-strain behaviors
include("test_plastic.jl")
include("test_viscoplastic.jl")
include("test_viscoelastic.jl")
include("test_differentiate.jl")
include("test_crystal_plasticity.jl")

# Test finite strain behaviors
include("test_hyperelastic.jl")
include("test_finitestrainplastic.jl")

# Test wrappers
include("test_rotated.jl")
