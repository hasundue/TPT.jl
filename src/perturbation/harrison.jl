"""
harrison.jl

Types and functions for perturbation using Wills-Harrison's tight-binding
approach.

References:
"""

immutable WHTB <: TBPerturbation
  c::Vector{Float64} # composition
  T::Float64 # temperature
  n::Float64 # coordination number
  zd::Vector{Float64} # number of d-electrons
  rd::Vector{Float64} # d-state radius
end

function WHTB(T::Number, zd::Number, rd::Number)
  WHTB([1.0], T, 12.0, [zd], [rd])
end


function pairpotential(whtb::WHTB)
  @attach(whtb, c, n, zd, rd)

  N = length(c)
  u = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    z̄d = (zd[i] + zd[j]) / 2

    u[i,j] = r -> -28.1/π * sqrt(12/n) * z̄d * (1 - z̄d/10) * (rd[i]*rd[j])^(3/2) / r^5 + 225/π^2 * z̄d * (rd[i]*rd[j])^3 / r^8
  end

  return u
end
