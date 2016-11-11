"""
nullpert.jl

Null Perturbation
"""

immutable NullPerturbation <: Perturbation
end

#
# Thermodynamic functions
#
kinetic(null::NullPerturbation, ref::ReferenceSystem)::Float64 = 0
entropy(null::NullPerturbation, ref::ReferenceSystem)::Float64 = 0
helmholtz(null::NullPerturbation, ref::ReferenceSystem)::Float64 = 0
