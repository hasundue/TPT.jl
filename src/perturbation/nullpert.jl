"""
nullpert.jl

Null Perturbation
"""

immutable NullPerturbation <: Perturbation
  ncomp::Int # number of components
end

#
# Thermodynamic functions
#
kinetic(null::NullPerturbation, ref::ReferenceSystem)::Float64 = 0
internal(null::NullPerturbation, ref::ReferenceSystem)::Float64 = 0
entropy(null::NullPerturbation, ref::ReferenceSystem, T::Float64)::Float64 = 0
helmholtz(null::NullPerturbation, ref::ReferenceSystem)::Float64 = 0

function pairpotential(null::NullPerturbation)::Array{Function,2}
  N::Int = null.ncomp
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue
    u(r)::Float64 = 0
    ret[i,j] = ret[j,i] = u
  end

  return ret
end
