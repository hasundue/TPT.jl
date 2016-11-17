"""
harrison.jl

Types and functions for perturbation using Wills-Harrison's tight-binding
approach.

References:
"""

immutable WHTB <: TBPerturbation
  zd::Vector{Float64} # number of d-electrons
  rd::Vector{Float64} # d-state radius
end

WHTB(zd::Real, rd::Real) = WHTB([zd], [rd])

ncomp(whtb::WHTB)::Int = length(whtb.zd)

function pairpotential(whtb::WHTB)::Array{Function,2}
  @attach(whtb, zd, rd)

  N = ncomp(whtb)
  ret = Array{Function,2}(N,N)
  n = 12 # coordination number

  for i in 1:N, j in 1:N
    i > j && continue

    z̄d = (zd[i] + zd[j]) / 2

    u(r) = -28.1/π * sqrt(12/n) * z̄d * (1 - z̄d/10) * (rd[i]*rd[j])^(3/2) / r^5 + 225/π^2 * z̄d * (rd[i]*rd[j])^3 / r^8

    ret[i,j] = ret[j,i] = u
  end

  return ret
end

function entropy(whtb::WHTB, ref::ReferenceSystem, T::Float64)
  @attach(whtb, zd, rd)

  N::Int = ncomp(ref)
  c::Vector{Float64} = composition(ref)
  d::Array{Float64,2} = nndistance(ref)

  S_el::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue

    Wd = 12 * 28.1/π * (rd[i]*rd[j])^(3/2) / d[i,j]^5
    nd = 10 / Wd

    a = i == j ? 1 : 2
    S_el += π^2 / 3 * kB^2 * T * (a * c[i]*c[j] * nd) / kB
  end

  return S_el
end

function internal(whtb::WHTB, ref::ReferenceSystem)::Float64
  U = 0
end
