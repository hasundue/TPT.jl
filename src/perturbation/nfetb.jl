"""
nfetb.jl

NFE-TB Perturbation
"""

immutable NFETB{Tn <: NFEPerturbation, Tt <: TBPerturbation} <: Perturbation
  nfe::Tn # NFE-like perturbation
  tb::Tt # TB-like perturbation
end

function pairpotential(nfetb::NFETB)::Array{Function,2}
  u_nfe = pairpotential(nfetb.nfe)
  u_tb = pairpotential(nfetb.tb)

  N = size(u_nfe, 1)
  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue
    u(r) =  u_nfe[i,j](r) + u_tb[i,j](r)
    ret[i,j] = ret[j,i] = u
  end

  return ret
end

function entropy(nfetb::NFETB, ref::ReferenceSystem, T::Float64)::Float64
  S_nfe = entropy(nfetb.nfe, ref, T)
  S_tb = entropy(nfetb.tb, ref, T)

  S = S_nfe + S_tb
end

function internal(nfetb::NFETB, ref::ReferenceSystem)::Float64
  U_nfe = internal(nfetb.nfe, ref)
  U_tb = internal(nfetb.tb, ref)

  U = U_nfe + U_tb
end
