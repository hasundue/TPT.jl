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
