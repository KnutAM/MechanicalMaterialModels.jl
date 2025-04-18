using Tensors, ForwardDiff, FiniteDiff

using MechanicalMaterialModels
import MechanicalMaterialModels as MechMat

using MaterialModelsBase
import MaterialModelsBase as MMB

using Test

include("utilities4testing.jl")

include("test_utils.jl")
include("test_hyperelastic.jl")

include("test_elastic.jl")
include("test_yield.jl")
include("test_hardening.jl")
include("test_overstress.jl")
include("test_plastic.jl")
include("test_viscoplastic.jl")
include("test_viscoelastic.jl")
include("test_differentiate.jl")

include("test_finitestrainplastic.jl")
