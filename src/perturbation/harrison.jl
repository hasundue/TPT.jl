"""
harrison.jl

Types and functions for perturbation using Wills-Harrison's tight-binding
approach.

References:
"""

immutable WHTB <: TBPerturbation
  zd::Vector{Float64} # number of d-electrons
  rd::Vector{Float64} # d-state radius
  c::Vector{Float64} # composition
end

WHTB(zd::Real, rd::Real) = WHTB([zd], [rd], [1.0])

ncomp(whtb::WHTB)::Int = length(whtb.zd)

function bandwidth(whtb::WHTB, ref::ReferenceSystem)::Float64
  @attach(whtb, zd, rd, c)

  N::Int = ncomp(ref)
  ρ::Float64 = totalnumberdensity(ref)
  g::Array{Function,2} = paircorrelation(ref)
  rcut::Array{Float64,2} = cutoffradius(ref)

  Wd²::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue
    n = i < j ? 1 : 2

    Vd(r) = 28.06/π * (rd[i]*rd[j])^(3/2) / r^5
    f(r) = Vd(r)^2 * g[i,j](r) * r^2
    integral = spline_integral(f, rcut[i,j], R_MAX, 64)

    Wd² += 12 * 4π*ρ * n*c[i]*c[j] * integral
  end

  Wd = √Wd²
end

function pairpotential(whtb::WHTB)::Array{Function,2}
  @attach(whtb, zd, rd, c)

  N::Int = ncomp(whtb)

  r̄d = prod(rd .^ c)
  γ = 12 # coordination number

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    z̄d = (zd[i] + zd[j]) / 2

    u(r) = -28.1/π * sqrt(12/γ) * z̄d * (1 - z̄d/10) * r̄d^3 / r^5 +
           225/π^2 * z̄d * (rd[i]*rd[j])^3 / r^8

    ret[i,j] = ret[j,i] = u
  end

  return ret
end

function entropy(whtb::WHTB, Wd::Float64, T::Float64)::Float64
  nd = 10 / Wd
  S_el = π^2 / 3 * kB^2 * T * nd / kB
end

function entropy(whtb::WHTB, ref::ReferenceSystem, T::Float64)::Float64
  entropy(whtb, bandwidth(whtb, ref), T)
end

function internal(whtb::WHTB, ref::ReferenceSystem)::Float64
  U = 0
end
