var documenterSearchIndex = {"docs":
[{"location":"utility_functions/","page":"Utility Functions","title":"Utility Functions","text":"CurrentModule = MechanicalMaterialModels","category":"page"},{"location":"utility_functions/#Utility-functions","page":"Utility Functions","title":"Utility functions","text":"","category":"section"},{"location":"utility_functions/","page":"Utility Functions","title":"Utility Functions","text":"In the course of defining the material behavior, a set of functions that may be useful outside the package internals are documented here, but these are not exported. ","category":"page"},{"location":"utility_functions/","page":"Utility Functions","title":"Utility Functions","text":"vonmises\nvonmises_and_gradient\nmacaulay\nconvert_hooke_param","category":"page"},{"location":"utility_functions/#MechanicalMaterialModels.vonmises","page":"Utility Functions","title":"MechanicalMaterialModels.vonmises","text":"function vonmises(σ::SecondOrderTensor{3})\n\nCalculate the von Mises effective stress for a 2nd order tensor\n\n\n\n\n\n","category":"function"},{"location":"utility_functions/#MechanicalMaterialModels.vonmises_and_gradient","page":"Utility Functions","title":"MechanicalMaterialModels.vonmises_and_gradient","text":"function vonmises_and_gradient(σ::SecondOrderTensor{3})\n\nCalculate the von Mises effective stress for a 2nd order tensor as well as the gradient, ν\n\n\n\n\n\n","category":"function"},{"location":"utility_functions/#MechanicalMaterialModels.macaulay","page":"Utility Functions","title":"MechanicalMaterialModels.macaulay","text":"function macaulay(x)\n\nCalculate the macaulay bracket of x, langle x rangle\n\nlangle x rangle = leftlbrace beginmatrix 0  xleq 0  x  x0 endmatrix right \n\n\n\n\n\n","category":"function"},{"location":"utility_functions/#MechanicalMaterialModels.convert_hooke_param","page":"Utility Functions","title":"MechanicalMaterialModels.convert_hooke_param","text":"convert_hooke_param(T::Symbol; p1, p2)\n\nConvert the hooke (isotropic) parameters p1 and p2 to the parameter with symbol T, e.g. G = convert_hooke_param(:G; E = 210e3, ν = 0.3)\n\n\n\n\n\n","category":"function"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"CurrentModule = MechanicalMaterialModels","category":"page"},{"location":"finite_strains/#Models-for-Finite-Strains","page":"Finite Strains","title":"Models for Finite Strains","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"Implementations of mechanical (stress-strain) material models following  the MaterialModelsBase.jl interface. ","category":"page"},{"location":"finite_strains/#Hyperelasticity","page":"Finite Strains","title":"Hyperelasticity","text":"","category":"section"},{"location":"finite_strains/#Neo-Hookean-Materials","page":"Finite Strains","title":"Neo-Hookean Materials","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"NeoHooke\nCompressibleNeoHooke","category":"page"},{"location":"finite_strains/#MechanicalMaterialModels.NeoHooke","page":"Finite Strains","title":"MechanicalMaterialModels.NeoHooke","text":"NeoHooke(; G)\n\nThe incompressible neo-Hookean formulation with shear modulus G defined by the potential\n\nvarPsi(boldsymbolC) = fracG2 left fracmathrmtr(boldsymbolC)sqrt3det(boldsymbolC) - 3right\n\nwhere boldsymbolC is the Right Cauchy-Green deformation tensor. \n\nNote that sometimes, the division by sqrt3det(boldsymbolC) is omitted,  which has no influence if C truly is incompressible, i.e. det(boldsymbolC) = 1. \n\n\n\n\n\n","category":"type"},{"location":"finite_strains/#MechanicalMaterialModels.CompressibleNeoHooke","page":"Finite Strains","title":"MechanicalMaterialModels.CompressibleNeoHooke","text":"CompressibleNeoHooke(; G, K)\n\nA compressible neo-Hookean formulation defined by the potential\n\nvarPsi(boldsymbolC) = \nfracG2 left fracmathrmtr(boldsymbolC)sqrt3det(boldsymbolC) - 3right\n+ varPsi_mathrmvolleft(sqrtdet(boldsymbolC)right) \nquad \nvarPsi_mathrmvol(J) = fracK2 left J - 1 right^2\n\nNote that there are many different variations of this model considering the volumetric part, varPsi_mathrmvol.\n\n\n\n\n\n","category":"type"},{"location":"finite_strains/#finite_strain_plasticity","page":"Finite Strains","title":"Plasticity","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"FiniteStrainPlastic","category":"page"},{"location":"finite_strains/#MechanicalMaterialModels.FiniteStrainPlastic","page":"Finite Strains","title":"MechanicalMaterialModels.FiniteStrainPlastic","text":"FiniteStrainPlastic(;elastic, yield, isotropic, kinematic, overstress)\n\nA finite-strain plasticity model with modular elastic laws, yield criteria,  multiple isotropic and kinematic hardening contributions, and either rate-independent or viscoplastic response.\n\nKeyword arguments\n\nelastic::AbstractHyperElastic\nHyperelastic law, see e.g. CompressibleNeoHooke\nyield::YieldCriterion\nYield criterion, including the initial yield limit. If yield::Real is given, VonMises(yield) is used. \nisotropic::Union{AbstractIsotropicHardening,Tuple}\nIsotropic hardening laws, see e.g. Voce\nkinematic::Union{AbstractKinematicHardening,Tuple}\nKinematic hardening laws, see e.g. ArmstrongFrederick\noverstress::Union{RateIndependent,OverstressFunction}\nRate dependence, see e.g. NortonOverstress\nDefaults to RateIndependent()\n\nExample\n\nm = FiniteStrainPlastic(elastic = CompressibleNeoHooke(G=80.e3, K=160.e3),\n            yield = 300.0,\n            isotropic = (Voce(Hiso=-100.e3, κ∞=-100.0),Voce(Hiso=10.e3, κ∞=200.0)),\n            kinematic = (ArmstrongFrederick(Hkin=200.e3, β∞=200.0),\n                         OhnoWang(Hkin=1000.e3, β∞=200.0, m=3.0)),\n            overstress = NortonOverstress(;tstar=1.0, nexp=2.0))\n\n\n\n\n\n","category":"type"},{"location":"finite_strains/#fsp_theory","page":"Finite Strains","title":"Theory","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"The overall modeling framework and general equations for FiniteStrainPlastic is described below. The exact model is specified by the provided laws as keyword arguments during construction.","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"The model is based on the multiplicative split of the deformation gradient, boldsymbolF,  into elastic, boldsymbolF_mathrme, and plastic, boldsymbolF_mathrmp, parts,","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"boldsymbolF=boldsymbolF_mathrme boldsymbolF_mathrmp","category":"page"},{"location":"finite_strains/#fsp_elasticity","page":"Finite Strains","title":"Elastic law","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"The 2nd Piola-Kirchhoff stress, boldsymbolS_mathrme,  on the intermediate (elastic) configuration is then  calculated from the hyperelastic potential, varPsi(boldsymbolC_mathrme),","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"boldsymbolS_mathrme = 2fracpartial varPsipartial boldsymbolC_mathrme","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"where the elastic Right Cauchy-Green deformation tensor, boldsymbolC_mathrme = boldsymbolF_mathrme^mathrmTboldsymbolF_mathrme. We further introduce the  Mandel stress as","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"boldsymbolM = 2 boldsymbolC_mathrme fracpartial varPsipartial boldsymbolC_mathrme","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"where boldsymbolM = boldsymbolF_mathrme^mathrmT boldsymboltau boldsymbolF_mathrme^mathrm-T and boldsymboltau is the Kirchhoff stress and boldsymbolP = boldsymbolF_mathrme boldsymbolS_mathrme boldsymbolF_mathrmp^-T is the 1st Piola-Kirchhoff stress.","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"In summary, the elastic response is determined by the potential, varPsi(boldsymbolC_mathrme),  specified by the elastic law.","category":"page"},{"location":"finite_strains/#fsp_yield","page":"Finite Strains","title":"Yield Criterion","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"A yield criterion of the type ","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"varPhi = fleft( boldsymbolM - boldsymbolM_mathrmk right) - leftY_0 + kapparight","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"is formulated in terms of Mandel stresses, where  boldsymbolM_mathrmk is the total back-stress, and kappa the total isotropic hardening stress. ","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"boldsymbolM_mathrmk = sum_i=1^N_mathrmkin boldsymbolM_mathrmki\nquad\nkappa = sum_i=1^N_mathrmiso kappa_i","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"How the back-stresses, boldsymbolM_mathrmki, and isotropic hardening stresses, kappa_i, are calculated is further described below.","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"The initial yield limit, Y_0, is passed to the yield criterion along with potentially other parameters affecting the function f. ","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"In this framework, an associative plastic flow is used to obtain the plastic deformations, specifically","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"boldsymbolL_mathrmp = dotboldsymbolF_mathrmp boldsymbolF_mathrmp^-1 = dotlambda mathrmnu quad mathrmnu = fracpartial varPhipartial boldsymbolM","category":"page"},{"location":"finite_strains/#fsp_isohard","page":"Finite Strains","title":"Isotropic Hardening","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"The isotropic hardening is formulated as kappa_i = g_mathrmisoi(lambda), determined by the specific isotropic hardening laws, see Isotropic hardening.","category":"page"},{"location":"finite_strains/#fsp_kinhard","page":"Finite Strains","title":"Kinematic Hardening","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"The back-stresses, boldsymbolM_mathrmki, are calculated as","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"boldsymbolM_mathrmki = 2boldsymbolc_mathrmki fracpartialvarPsi_mathrmkipartial boldsymbolc_mathrmki","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"where the deformation tensor associated with kinematic hardening,  boldsymbolc_mathrmki = boldsymbolF_mathrmki^-mathrmTboldsymbolF_mathrmki^-1.  The free energy for the kinematic hardening stress, varPsi_mathrmki, is always given by  the incompressible NeoHooke formulation, with modulus G=H3, where H is the hardening modulus given to the kinematic hardening law. The factor 13 gives the same initial plastic stiffness for both isotropic and kinematic hardening.","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"The evolution of the kinematic hardening deformation gradient, boldsymbolF_mathrmki, is then given by the specific kinematic hardening law,","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"boldsymbolL_mathrmp = dotboldsymbolF_mathrmp boldsymbolF_mathrmp^-1 = dotlambda boldsymbolg_mathrmki(boldsymbolnu boldsymbolM_mathrmki)","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"where boldsymbolg_mathrmki are determined by the specific  hardening laws, see Kinematic hardening.","category":"page"},{"location":"finite_strains/#fsp_ratedep","page":"Finite Strains","title":"Rate Dependence","text":"","category":"section"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"If overstress=RateIndependent(), the plastic multiplier, lambda, is obtained via the KKT-conditions,","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"dotlambda geq 0 quad varPhi leq 0 quad dotlambdavarPhi = 0","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"Otherwise, the overstress function, eta(varPhi), determines the evolution of lambda as ","category":"page"},{"location":"finite_strains/","page":"Finite Strains","title":"Finite Strains","text":"dotlambda = eta(varPhi Y_0 + kappa)","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"CurrentModule = MechanicalMaterialModels","category":"page"},{"location":"model_components/#Model-Components","page":"Model Components","title":"Model Components","text":"","category":"section"},{"location":"model_components/#Yield-criteria","page":"Model Components","title":"Yield criteria","text":"","category":"section"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"VonMises\nDruckerPrager","category":"page"},{"location":"model_components/#MechanicalMaterialModels.VonMises","page":"Model Components","title":"MechanicalMaterialModels.VonMises","text":"VonMises(; Y0)\n\nCreate a von Mises yield criterion with initial yield limit, Y_0, as Y0. The yield criterion is then defined as\n\nPhi = sqrtfrac32 left textdev left( boldsymbolsigma_mathrmred right) right - left Y_0 + Delta Y right = 0\n\nwhere boldsymbolsigma_mathrmred is the reduced (by kinematic hardening) stress tensor, and Delta Y the change of the initial  yield limit due to isotropic hardening (i.e. kappa).\n\n\n\n\n\n","category":"type"},{"location":"model_components/#MechanicalMaterialModels.DruckerPrager","page":"Model Components","title":"MechanicalMaterialModels.DruckerPrager","text":"DruckerPrager(; Y0, B)\n\nCreate a Drucker-Prager yield criterion, with initial yield limit, Y_0, as Y0, and pressure sensitivity B. The yield criterion is defined as \n\nPhi = sqrtfrac32 left mathrmdev left( boldsymbolsigma_mathrmred right) right - Bmathrmtrleft( boldsymbolsigma_mathrmred right) - left Y_0 + Delta Y right = 0\n\nwhere boldsymbolsigma_mathrmred is the reduced (by kinematic hardening) stress tensor, and Delta Y the change of the initial  yield limit due to isotropic hardening (i.e. kappa).\n\n\n\n\n\n","category":"type"},{"location":"model_components/#Isotropic-hardening","page":"Model Components","title":"Isotropic hardening","text":"","category":"section"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"Voce\nSwift","category":"page"},{"location":"model_components/#MechanicalMaterialModels.Voce","page":"Model Components","title":"MechanicalMaterialModels.Voce","text":"Voce(;Hiso, κ∞)\n\nExponentially saturating isotropic hardening\n\nkappa_i = g_mathrmisoi(lambda) = kappa_infty left1 - mathrmexpleft(fracH_mathrmisokappa_infty lambda right)right\n\nor alternatively as differential equations\n\ndotkappa_i = dotlambda H_mathrmiso left1 - frackappa_ikappa_inftyright\n\nArguments\n\nHiso: Isotropic hardening modulus, H_mathrmiso\nκ∞: Saturation hardening value, kappa_infty\n\n\n\n\n\n","category":"type"},{"location":"model_components/#MechanicalMaterialModels.Swift","page":"Model Components","title":"MechanicalMaterialModels.Swift","text":"Swift(; K, λ0, n)\n\nIsotropic hardening by the Swift power law\n\nkappa_i = g_mathrmisoi(lambda) = K leftlambda_0 + lambda right^n\n\nArguments\n\nK: K\nλ0: lambda_0\nn: n\n\n\n\n\n\n","category":"type"},{"location":"model_components/#Kinematic-hardening","page":"Model Components","title":"Kinematic hardening","text":"","category":"section"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"ArmstrongFrederick","category":"page"},{"location":"model_components/#MechanicalMaterialModels.ArmstrongFrederick","page":"Model Components","title":"MechanicalMaterialModels.ArmstrongFrederick","text":"ArmstrongFrederick(; Hkin, β∞)\n\nArmstrong-Frederick kinematic hardening law with modulus Hkin and saturation stress β∞. \n\n\n\n\n\n","category":"type"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"The evolution equation according to Armstrong-Frederick hardening  for boldsymbolb is given as","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"g_mathrmkini(nu boldsymbolbeta_i) = -boldsymbolnu + frac3boldsymbolbeta_i^mathrmT2beta_infty","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"noting that for small strains, boldsymbolbeta_i = boldsymbolbeta_i^mathrmT is symmetric. ","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"Temp ref: doi: 10.1179/096034007X207589","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"Delobelle","category":"page"},{"location":"model_components/#MechanicalMaterialModels.Delobelle","page":"Model Components","title":"MechanicalMaterialModels.Delobelle","text":"Delobelle(; Hkin, β∞, δ)\n\nKinematic hardening law according to Delobelle with hardening modulus Hkin, saturation stress, β∞, and scaling parameter δ, which scales between pure Armstrong-Frederick hardening, δ=1, and Burlet-Cailletaud hardening, δ=0.\n\n\n\n\n\n","category":"type"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"The evolution equation according to Delobelle hardening  for boldsymbolb is given as","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"g_mathrmkini(nu boldsymbolbeta_i) \n= - boldsymbolnu + deltafrac3boldsymbolbeta_i^mathrmT2beta_infty\n+ left1 - deltaright fracboldsymbolnuboldsymbolbeta_ibeta_inftyboldsymbolnu","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"noting that for small strains, boldsymbolbeta_i = boldsymbolbeta_i^mathrmT is symmetric. ","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"Temp ref: doi: 10.1016/S0749-6419(95)00001-1","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"OhnoWang","category":"page"},{"location":"model_components/#MechanicalMaterialModels.OhnoWang","page":"Model Components","title":"MechanicalMaterialModels.OhnoWang","text":"OhnoWang(; Hkin, β∞, m)\n\nKinematic hardening law according to Ohno-Wang with hardening Hkin, saturation stress, β∞, and exponent, m.\n\n\n\n\n\n","category":"type"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"The evolution equation according to Ohno-Wang hardening for boldsymbolb is given as","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"g_mathrmkini(nu boldsymbolbeta_i) \n= - boldsymbolnu \n  + frac3boldsymbolbeta_i^mathrmT2beta_infty \n    fraclangle boldsymbolnuboldsymbolbeta_i ranglebeta_infty\n    leftfracbeta_i^mathrmvMbeta_inftyright^m","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"where langle x rangle = mathrmmax(x 0) is the Macaulay bracket and beta_i^mathrmvM is the effective von Mises backstress, see vonmises. Note that for small strains, boldsymbolbeta_i = boldsymbolbeta_i^mathrmT is symmetric. ","category":"page"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"Temp ref: 10.1016/0749-6419(93)90042-O","category":"page"},{"location":"model_components/#Overstress-functions","page":"Model Components","title":"Overstress functions","text":"","category":"section"},{"location":"model_components/","page":"Model Components","title":"Model Components","text":"RateIndependent\nNortonOverstress","category":"page"},{"location":"model_components/#MechanicalMaterialModels.RateIndependent","page":"Model Components","title":"MechanicalMaterialModels.RateIndependent","text":"RateIndependent()\n\nThe evolution of the plastic multiplier for a rate-dependent material is given by  the so-called KKT loading/unloading conditions\n\ndotlambda geq 0 quad varPhi leq 0 quad dotlambdavarPhi = 0\n\n\n\n\n\n","category":"type"},{"location":"model_components/#MechanicalMaterialModels.NortonOverstress","page":"Model Components","title":"MechanicalMaterialModels.NortonOverstress","text":"NortonOverstress(; tstar, nexp)\n\nThe norton overstress function is defined as \n\neta(varPhi sigma_mathrmy) = frac1t_* leftlangle fracvarPhisigma_mathrmy rightrangle^n\n\nwhere the material parameters t_* (tstar) and n (nexp) represent the  relaxation time and overstress sensitivty.  \n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MechanicalMaterialModels","category":"page"},{"location":"#MechanicalMaterialModels","page":"Home","title":"MechanicalMaterialModels","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Implementations of mechanical (stress-strain) material models following  the MaterialModelsBase.jl  interface, and to use this package that and  Newton.jl must be installed:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(url=\"https://github.com/KnutAM/MaterialModelsBase.jl\")\nPkg.add(url=\"https://github.com/KnutAM/Newton.jl\")\nPkg.add(url=\"https://github.com/KnutAM/MechanicalMaterialModels.jl\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"See MaterialModelsBase.jl's documentation for details on how to use materials.","category":"page"},{"location":"#Examples-of-available-material-models","page":"Home","title":"Examples of available material models","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"LinearElastic: Isotropic and anisotropic linear elasticity for small strains\nGeneralizedMaxwell: The standard generalized linear viscoelasticity model with arbitrarily many maxwell-chains.\nPlastic: Small strain plasticity supporting\nDifferent yield surfaces, e.g. von Mises and Drucker-Prager\nDifferent isotropic and kinematic hardening laws\nBoth viscoplasticity with overstress functions and standard rate-indepenent plasticity.\nCompressibleNeoHooke: One of several hyperelasticity models included for finite strains\nFiniteStrainPlastic: Finite strain plasticity with the same features as Plastic, but for finite strains.","category":"page"},{"location":"devdocs/","page":"Developer Docs","title":"Developer Docs","text":"CurrentModule = MechanicalMaterialModels","category":"page"},{"location":"devdocs/#Developer-Documentation","page":"Developer Docs","title":"Developer Documentation","text":"","category":"section"},{"location":"devdocs/","page":"Developer Docs","title":"Developer Docs","text":"The following functions and their docstrings may be helpful  for understanding the code, but these may change.","category":"page"},{"location":"devdocs/","page":"Developer Docs","title":"Developer Docs","text":"maketuple_or_nothing\nget_promoted_type\nbaseof\nvector_residual!\ncompute_potential\ncompute_stress\ncompute_tangent\ncompute_stress_and_tangent\nstatic_vector\nget_resid_eltype\nget_num_unknowns","category":"page"},{"location":"devdocs/#MechanicalMaterialModels.maketuple_or_nothing","page":"Developer Docs","title":"MechanicalMaterialModels.maketuple_or_nothing","text":"maketuple_or_nothing(x)\n\nx is a single value: Convert to a Tuple of length 1\nx is a Tuple or Nothing: Return x\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.get_promoted_type","page":"Developer Docs","title":"MechanicalMaterialModels.get_promoted_type","text":"get_promoted_type(args...)\n\nGet the promoted type for the type of the arguments,  e.g. get_promoted_type(1, 1.f0) is Float32\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.baseof","page":"Developer Docs","title":"MechanicalMaterialModels.baseof","text":"baseof(t::AbstractTensor)\n\nGet the base-type of t, i.e. if t::SymmetricTensor{2,3,Float64,6}, baseof(t) returns SymmetricTensor{2,3}\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.vector_residual!","page":"Developer Docs","title":"MechanicalMaterialModels.vector_residual!","text":"vector_residual!(rf::Function, r_vector::AbstractVector, x_vector::AbstractVector, x::AbstractResidual)\n\nMakes it easy to construct a mutating vector residual function from a tensor-like equation, r = rf(x) = residual(x, args...), e.g. rf!(r_vector, x_vector) = vector_residual!(z -> residual(z, args...), r_vector, x_vector, x)\n\nThe input x and output r of rf should have the same type, RT <: AbstractResidual,  and support MaterialModelsBase.tovector! and MaterialModelsBase.fromvector!\n\nThe approach was adopted from https://github.com/kimauth/MaterialModels.jl\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.compute_potential","page":"Developer Docs","title":"MechanicalMaterialModels.compute_potential","text":"compute_potential(m::AbstractHyperElastic, C::SymmetricTensor)\n\nCompute the potential Ψ(C) given the Right-Cauchy-Green deformation tensor C  for the hyperelastic material m.\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.compute_stress","page":"Developer Docs","title":"MechanicalMaterialModels.compute_stress","text":"compute_stress(m::AbstractHyperElastic, C::SymmetricTensor)\n\nCompute the 2nd PiolaKirchhoff stress, S = 2 ∂Ψ/∂C, for the potential Ψ defined by m for the Right-Cauchy-Green deformation tensor C \n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.compute_tangent","page":"Developer Docs","title":"MechanicalMaterialModels.compute_tangent","text":"compute_tangent(m::AbstractHyperElastic, C::SymmetricTensor)\n\nCompute the tangent stiffness, ∂S/∂E = 4 ∂²Ψ/∂C², for the potential Ψ defined by m for the Right-Cauchy-Green deformation tensor C.\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.compute_stress_and_tangent","page":"Developer Docs","title":"MechanicalMaterialModels.compute_stress_and_tangent","text":"compute_stress_and_tangent(m::AbstractHyperElastic, C::SymmetricTensor)\n\nCompute both the 2nd Piola-Kirchhoff stress, S, and the tangent stiffness, ∂S/∂E = 4 ∂²Ψ/∂C²,  for the potential Ψ defined by m for the Right-Cauchy-Green deformation tensor C. This is normally more efficient than computing the stress and tangents invidividually.\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.static_vector","page":"Developer Docs","title":"MechanicalMaterialModels.static_vector","text":"static_vector(args::Union{Number, SVector}...)\n\nConvert the elements of args into an SVector.\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.get_resid_eltype","page":"Developer Docs","title":"MechanicalMaterialModels.get_resid_eltype","text":"get_resid_eltype(res::AbstractResidual)\n\nReturn the element type used to store res as a vector\n\n\n\n\n\n","category":"function"},{"location":"devdocs/#MechanicalMaterialModels.get_num_unknowns","page":"Developer Docs","title":"MechanicalMaterialModels.get_num_unknowns","text":"get_num_unknowns(res::AbstractResidual)\n\nReturn the number of unknowns for the residual res\n\n\n\n\n\n","category":"function"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"CurrentModule = MechanicalMaterialModels","category":"page"},{"location":"small_strains/#Models-for-Small-Strains","page":"Small Strains","title":"Models for Small Strains","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"Implementations of mechanical (stress-strain) material models following  the MaterialModelsBase.jl interface. ","category":"page"},{"location":"small_strains/#Elasticity","page":"Small Strains","title":"Elasticity","text":"","category":"section"},{"location":"small_strains/#Linear-Isotropic-Elasticity","page":"Small Strains","title":"Linear Isotropic Elasticity","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"LinearElastic(::Val{:isotropic})","category":"page"},{"location":"small_strains/#MechanicalMaterialModels.LinearElastic-Tuple{Val{:isotropic}}","page":"Small Strains","title":"MechanicalMaterialModels.LinearElastic","text":"LinearElastic(; E, ν)\nLinearElastic{:isotropic}(; E, ν)\n\nCreate an isotropic LinearElastic material with Young's modulus, E, and Poisson's ratio ν, such that\n\nboldsymbolsigma = 2mu boldsymbolepsilon + lambda mathrmtr(boldsymbolepsilon) boldsymbolI \n\nwhere the Lamé parameters, mu and lambda are defined as\n\nmu = fracE2(1+nu) quad lambda=fracEnu(1+nu)(1-2nu)\n\n\n\n\n\n","category":"method"},{"location":"small_strains/#Linear-Anisotropic-Elasticity","page":"Small Strains","title":"Linear Anisotropic Elasticity","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"LinearElastic(::Val{:general})\nLinearElastic(::Val{:cubicsymmetry})","category":"page"},{"location":"small_strains/#MechanicalMaterialModels.LinearElastic-Tuple{Val{:general}}","page":"Small Strains","title":"MechanicalMaterialModels.LinearElastic","text":"LinearElastic(C::SymmetricTensor{4,3})\nLinearElastic{:general}(C::SymmetricTensor{4,3})\n\nCreate a general LinearElastic material with the 4th order elastic stiffness tensor boldsymbolC, such that  boldsymbolsigma = boldsymbolCboldsymbolepsilon. \n\n\n\n\n\n","category":"method"},{"location":"small_strains/#MechanicalMaterialModels.LinearElastic-Tuple{Val{:cubicsymmetry}}","page":"Small Strains","title":"MechanicalMaterialModels.LinearElastic","text":"LinearElastic{:cubicsymmetry}(; C1111::T, C1122::T, C1212::T) where {T}\n\nCreate a LinearElastic material where the stiffness tensor, boldsymbolC, possesses cubic symmetry along the coordinate axes. Using the 9-component Voigt notation, boldsymbolC can be expressed as\n\nboldsymbolC = \nbeginbmatrix\nC_1111  C_1122  C_1122  0  0  0  0  0  0 \nC_1122  C_1111  C_1122  0  0  0  0  0  0 \nC_1122  C_1122  C_1111  0  0  0  0  0  0 \n0  0  0  2C_1212  0  0  0  0  0 \n0  0  0  0  2C_1212  0  0  0  0 \n0  0  0  0  0  2C_1212  0  0  0 \n0  0  0  0  0  0  2C_1212  0  0 \n0  0  0  0  0  0  0  2C_1212  0 \n0  0  0  0  0  0  0  0  2C_1212 \nendbmatrix\n\n\n\n\n\n","category":"method"},{"location":"small_strains/#small_strain_viscoelasticity","page":"Small Strains","title":"Viscoelasticity","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"GeneralizedMaxwell\nMaxwell","category":"page"},{"location":"small_strains/#MechanicalMaterialModels.GeneralizedMaxwell","page":"Small Strains","title":"MechanicalMaterialModels.GeneralizedMaxwell","text":"GeneralizedMaxwell(elastic::LinearElastic, chains::Maxwell...)\n\nCreate a generalized Maxwell model with an arbitrary number of Maxwell chains. The elastic part refers to the long-term stiffness contribution. \n\nCurrently, GeneralizedMaxwell and Maxwell are specific to isotropic viscous behavior, this should be generalized.\n\n\n\n\n\n","category":"type"},{"location":"small_strains/#MechanicalMaterialModels.Maxwell","page":"Small Strains","title":"MechanicalMaterialModels.Maxwell","text":"Maxwell(; G, η, t)\n\nCreate a maxwell chain element with shear stiffness G. Either the viscosity, η,  or the relaxation time, t = η / G, should be supplied (not both).\n\n\n\n\n\n","category":"type"},{"location":"small_strains/#small_strain_plasticity","page":"Small Strains","title":"Plasticity","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"Plastic","category":"page"},{"location":"small_strains/#MechanicalMaterialModels.Plastic","page":"Small Strains","title":"MechanicalMaterialModels.Plastic","text":"Plastic(;elastic, yield, isotropic, kinematic, overstress)\n\nA small-strain plasticity model with modular elastic laws, yield criteria,  multiple isotropic and kinematic hardening contributions, and either rate-independent or viscoplastic response.\n\nKeyword arguments\n\nelastic::AbstractMaterial\nElastic law, see e.g. LinearElastic\nyield::YieldCriterion\nYield criterion, including the initial yield limit. If yield::Real is given, VonMises(yield) is used. \nisotropic::Union{AbstractIsotropicHardening,Tuple}\nIsotropic hardening laws, see e.g. Voce\nkinematic::Union{AbstractKinematicHardening,Tuple}\nKinematic hardening laws, see e.g. ArmstrongFrederick\noverstress::Union{RateIndependent,OverstressFunction}\nRate dependence, see e.g. NortonOverstress\nDefaults to RateIndependent()\n\nExample\n\nm = Plastic(elastic = LinearElastic(E=210.e3, ν=0.3),\n            yield = 100.0,\n            isotropic = (Voce(Hiso=-100.e3, κ∞=-100.0),Voce(Hiso=10.e3, κ∞=200.0)),\n            kinematic = (ArmstrongFrederick(Hkin=200.e3, β∞=200.0),\n                         OhnoWang(Hkin=1000.e3, β∞=200.0, m=3.0)),\n            overstress = NortonOverstress(;tstar=1.0, nexp=2.0))\n\n\n\n\n\n","category":"type"},{"location":"small_strains/#ssp_theory","page":"Small Strains","title":"Theory","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"While the exact model response is given by the laws in elastic, yield, isotropic, kinematic, and overstress, the generic model equations are described below. ","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"The stress is calculated from the elastic strains, boldsymbolepsilon_mathrme, obtained via the  additive decomposition, boldsymbolepsilon = boldsymbolepsilon_mathrme + boldsymbolepsilon_mathrmp.  The elastic law is specified by m.elastic and is evaluated by giving it the elastic strain. ","category":"page"},{"location":"small_strains/#ssp_yield","page":"Small Strains","title":"Yield Criterion","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"A yield criterion of the type ","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"varPhi = fleft( boldsymbolsigma - boldsymbolbeta right) - leftY_0 + kapparight","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"is assumed. Here, boldsymbolbeta = sum_i=1^N_mathrmkin boldsymbolbeta_i is the total back-stress,  and kappa = sum_i=1^N_mathrmiso kappa_i is the total isotropic hardening stress. The initial yield limit  is passed to the yield criterion along with potentially other parameters.  The evolution laws for boldsymbolbeta_i and kappa_i are given by the kinematic and isotropic hardening laws.","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"Associative plastic flow is used to obtain the plastic strains,","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"dotepsilon_mathrmp = dotlambda leftfracpartial fpartial boldsymbolsigmarightvert_left( boldsymbolsigma - boldsymbolbeta right)\n= dotlambda boldsymbolnu","category":"page"},{"location":"small_strains/#ssp_isohard","page":"Small Strains","title":"Isotropic Hardening","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"The isotropic hardening stress is kappa_i = -H_mathrmisoi k_i, where the evolution law for k_i is","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"dotk_i = dotlambda g_mathrmisoi(kappa_i)","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"where g_mathrmisoi is specified by an isotropic hardening law, see Isotropic hardening. The hardening modulus, H_mathrmisoi, is a parameter given to the isotropic hardening law. ","category":"page"},{"location":"small_strains/#ssp_kinhard","page":"Small Strains","title":"Kinematic Hardening","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"The back-stress is boldsymbolbeta_i = - 23H_mathrmkini boldsymbolb_i, where the evolution law for boldsymbolb_i is","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"dotboldsymbolb_i = dotlambda boldsymbolg_mathrmkini(boldsymbolnu boldsymbolbeta_i)","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"where g_mathrmkini is specified by a kinematic hardening law, see  Kinematic hardening. The hardening modulus, H_mathrmkini, is a parameter given to the kinematic hardening law. ","category":"page"},{"location":"small_strains/#ssp_ratedep","page":"Small Strains","title":"Rate Dependence","text":"","category":"section"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"If overstress=RateIndependent(), the plastic multiplier, lambda, is obtained via the KKT-conditions,","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"dotlambda geq 0 quad varPhi leq 0 quad dotlambdavarPhi = 0","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"Otherwise, the overstress function, eta(varPhi), determines the evolution of lambda as ","category":"page"},{"location":"small_strains/","page":"Small Strains","title":"Small Strains","text":"dotlambda = eta(varPhi (Y_0 + kappa))","category":"page"}]
}
