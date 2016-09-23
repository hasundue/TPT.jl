"""
harrison.jl

Types and functions for perturbation using Wills-Harrison's tight-binding
approach.

References:
"""

immutable WHTB <: Perturbation
  ρ::Float64 # number density
  c::Vector{Float64} # composition
  T::Float64 # temperature
  n::Float64 # coordination number
  zd::Vector{Float64} # number of d-electrons
  rd::Vector{Float64} # d-state radius
end

function WHTB(ρ::Number, T::Number, zd::Number, rd::Number)
  WHTB(ρ, [1.0], T, 12.0, [zd], [rd])
end


function pairpotential(whtb::WHTB)
  @attach(whtb, c, n, zd, rd)

  z̄d = sum(c .* zd)

  N = length(c)
  ub = Array{Function}(N,N)
  uc = Array{Function}(N,N)

  for k in 1:length(ub)
    i, j = ind2sub((N,N), k)

    ub[i,j] = r -> -28.1/π * sqrt(12/n) * z̄d * (1 - z̄d/10) * (rd[i]*rd[j])^(3/2) / r^5
    uc[i,j] = r -> 225/π^2 * z̄d * (rd[i]*rd[j])^3 / r^8
  end

  return ub, uc
end
