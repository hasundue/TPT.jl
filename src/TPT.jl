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

       # Thermodynamic properties
       kinetic,
       entropy,
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

       # utils
       spline,

       # constants
       kB

import Optim
import NLopt
using Dierckx

include("utils.jl")
include("types.jl")
include("tptsystem.jl")
include("constants.jl")
include(joinpath("reference", "reference.jl"))
include(joinpath("reference", "ahs.jl"))
include(joinpath("reference", "wca.jl"))
include(joinpath("perturbation", "nullpert.jl"))
include(joinpath("perturbation", "lennardjones.jl"))
include(joinpath("perturbation", "nfetb.jl"))
include(joinpath("perturbation", "nfe.jl"))
include(joinpath("perturbation", "ashcroft.jl"))
include(joinpath("perturbation", "bretonnet_silbert.jl"))
include(joinpath("perturbation", "harrison.jl"))

end # module
