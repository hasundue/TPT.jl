"""
harrison.jl

Types and functions for perturbation using Wills-Harrison's tight-binding
approach.

References:
"""

const η_original = (-16.2, 8.75, -2.39)
const η_modified = (-21.22, 12.60, -2.29)
const η_dict = Dict(:original => η_original, :modified => η_modified)

immutable WHTB <: TBPerturbation
  version::Symbol # :original or :modified parameters
  zd::Vector{Float64} # number of d-electrons
  rd::Vector{Float64} # d-state radius
  c::Vector{Float64} # composition
end

function WHTB(zd::Real, rd::Real; version=:original)
  WHTB(version, [zd], [rd], [1.0])
end

function WHTB(zd::Vector, rd::Vector, c::Vector; version=:original)
  WHTB(version, zd, rd, c)
end

ncomp(whtb::WHTB)::Int = length(whtb.zd)

function coreradius(whtb::WHTB)::Array{Float64,2}
  N::Int = ncomp(whtb)
  zeros(N,N)
end

function cutoffradius(whtb::WHTB)::Array{Float64,2}
  rd::Vector{Float64} = whtb.rd
  N::Int = ncomp(whtb)

  [ 10 * (rd[i] + rd[j]) / 2 for i in 1:N, j in 1:N ]
end

function transfercoefficient(whtb::WHTB)::Float64
  (η_ddσ, η_ddπ, η_ddδ) =
    get(η_dict, whtb.version, zeros(3))

  ( (η_ddσ^2 + 2η_ddπ^2 + 2η_ddδ^2) / 5 )^(1/2)
end

function transfermatrixelement(whtb::WHTB)::Array{Function,2}
  N::Int = ncomp(whtb)
  rd::Vector{Float64} = whtb.rd
  A::Float64 = transfercoefficient(whtb)

  [ V(r::Float64)::Float64 = A * (rd[i]*rd[j])^(3/2) / r^5
    for i in 1:N, j in 1:N ]
end

function overlapcoefficient(whtb::WHTB)::Float64
  (η_ddσ, η_ddπ, η_ddδ) =
    get(η_dict, whtb.version, zeros(3))

  (σ_ddσ, σ_ddπ, σ_ddδ) = (5/π, -5/π, 5/2π)

  -2/5 * (σ_ddσ*η_ddσ + 2σ_ddπ*η_ddπ + 2σ_ddδ*η_ddδ)
end

function overlapmatrixelement(whtb::WHTB)::Array{Function,2}
  N::Int = ncomp(whtb)
  rd::Vector{Float64} = whtb.rd
  B::Float64 = overlapcoefficient(whtb)

  [ SV(r::Float64)::Float64 = B * (rd[i]*rd[j])^3 / r^8
    for i in 1:N, j in 1:N ]
end

function bandwidth(whtb::WHTB, ref::ReferenceSystem)::Float64
  @attach(whtb, zd, rd, c)

  N::Int = ncomp(ref)
  ρ::Float64 = totalnumberdensity(ref)
  g::Array{Function,2} = paircorrelation(ref)
  rmin::Array{Float64,2} = emptyradius(ref)
  rcut::Array{Float64,2} = cutoffradius(whtb)
  V::Array{Function,2} = transfermatrixelement(whtb)

  Wd²::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue
    n = i == j ? 1 : 2

    gs = spline(g[i,j], rmin[i,j], rcut[i,j], 32)
    f(r) = V[i,j](r)^2 * gs(r) * r^2

    Wd² += 12 * 4π*ρ * n*c[i]*c[j] * ∫(f, rmin[i,j], rcut[i,j])
  end

  Wd = √Wd²
end

function pairpotential(whtb::WHTB)::Array{Function,2}
  @attach(whtb, zd, rd, c)

  N::Int = ncomp(whtb)
  A::Float64 = transfercoefficient(whtb)
  B::Float64 = overlapcoefficient(whtb)

  γ = 12 # coordination number

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    z̄d = (zd[i] + zd[j]) / 2
    r̄d = √(rd[i]*rd[j])

    u(r) = -A * sqrt(12/γ) * z̄d * (1 - z̄d/10) * r̄d^3 / r^5 +
           B * z̄d * r̄d^6 / r^8

    ret[i,j] = ret[j,i] = u
  end

  return ret
end

function pairpotential_derivative(whtb::WHTB)::Array{Function,2}
  @attach(whtb, zd, rd, c)

  N::Int = ncomp(whtb)
  A::Float64 = transfercoefficient(whtb)
  B::Float64 = overlapcoefficient(whtb)

  γ = 12 # coordination number

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    z̄d = (zd[i] + zd[j]) / 2
    r̄d = √(rd[i]*rd[j])

    u(r) = 5A * sqrt(12/γ) * z̄d * (1 - z̄d/10) * r̄d^3 / r^6 +
           -8B * z̄d * r̄d^6 / r^9

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
