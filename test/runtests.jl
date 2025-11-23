include("setup_tests.jl")

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
