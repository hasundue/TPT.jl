"""
nfetb.jl

NFE-TB Perturbation
"""

immutable NFETB{Tn <: NFEPerturbation, Tt <: TBPerturbation} <: Perturbation
  nfe::Tn # NFE-like perturbation
  tb::Tt # TB-like perturbation
end

function pairpotential(nfetb::NFETB)
  u_nfe = pairpotential(nfetb.nfe)
  u_tb = pairpotential(nfetb.tb)

  N = size(u_nfe, 1)
  u = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    u[i,j] = r -> u_nfe[i,j](r) + u_tb[i,j](r)
  end

  return u
end

function potentialenergy(ref::ReferenceSystem, nfetb::NFETB)
  g = prdf(ref)
  u = pairpotential(nfetb)
  ρ = numberdensity(ref)
  c = composition(ref)
  N = length(c)
  I = Array{Float64}(N,N)

  for i in 1:N, j in 1:N
    us = spline(u[i,j], R_MIN, R_MAX, 64)

    f = r -> g[i,j](r) == 0. ? 0. : g[i,j](r) * us(r) * r^2
    I[i,j] = c[i]*c[j] * ∫(f, R_MIN, R_MAX)
  end

  return 2π*ρ*sum(I)
end
