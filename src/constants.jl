"""
constants.jl

Physical constants and so on
"""

# Boltzmann constant
const kB = 3.167e-6 # a.u.

# cutoff parameters
const R_MIN = eps(Float64)
const R_MAX = 30.

const Q_MAX = 10.

# other internal constants
const InvalTemp = -1.
const InvalMass = Float64[]
