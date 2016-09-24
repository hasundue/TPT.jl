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
  u_b, u_c = pairpotential(nfetb.tb)

  N, M = size(u_nfe)
  u = Array{Function}(N,M)

  for i in 1:N, j in 1:M
    u[i,j] = r -> u_nfe[i,j](r) + u_b[i,j](r) + u_c[i,j](r)
  end

  return u
end
