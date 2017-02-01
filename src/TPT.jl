module TPT

export TPTSystem,

       # Basic information
       ncomp,
       composition,
       numberdensity,
       totalnumberdensity,
       temperature,

       # Structural properties
       structurefactor,
       paircorrelation,
       cavityfunction,
       nndistance,

       # Interatomic interaction
       pairpotential,
       pairpotential_minimum,
       pairpotential_minimizer,
       hsdiameter_estimate,

       # Thermodynamic properties
       kinetic,
       entropy,
       entropy_gas,
       entropy_conf,
       internal,
       helmholtz,

       # Reference systems
       AHS,
       WCA,

       # Hard-sphere system
       hsdiameter,
       packingfraction,
       totalpackingfraction,
       contactvalue,
       contactgradient,

       # Perturbation
       LennardJones,
       NFE,
       WHTB,
       NFETB,
       BOTB,

       # WCA
       blipfunction,

       # NFE
       Ashcroft,
       BretonnetSilbert,
       fermiwavenumber,
       formfactor,
       screenedformfactor,
       dielectric,
       localfiled,
       wnechar,

       # TB
       bandwidth,

       # utils
       spline,

       # constants
       kB

import Optim
import NLopt
using Dierckx
using Polynomials

include("utils.jl")
include("types.jl")
include("tptsystem.jl")
include("constants.jl")
include(joinpath("reference", "reference.jl"))
include(joinpath("reference", "abstractwca.jl"))
include(joinpath("reference", "ahs.jl"))
include(joinpath("reference", "wca.jl"))
include(joinpath("reference", "lwca.jl"))
include(joinpath("perturbation", "perturbation.jl"))
include(joinpath("perturbation", "nullpert.jl"))
include(joinpath("perturbation", "lennardjones.jl"))
include(joinpath("perturbation", "nfetb.jl"))
include(joinpath("perturbation", "nfe.jl"))
include(joinpath("perturbation", "ashcroft.jl"))
include(joinpath("perturbation", "bretonnet_silbert.jl"))
include(joinpath("perturbation", "harrison.jl"))
include(joinpath("perturbation", "hausleitner.jl"))

end # module
