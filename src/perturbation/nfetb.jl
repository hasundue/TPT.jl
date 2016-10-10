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

function internalenergy(ref::ReferenceSystem, nfetb::NFETB)
  g = prdf(ref)
  u = pairpotential(nfetb)
  ρ = numberdensity(ref)
  c = composition(ref)
  N = length(c)
  I = Array{Float64}(N,N)

  for i in 1:N, j in 1:N
    if i > j
      continue
    end

    I[i,j] = c[i]*c[j] * ∫(r -> u[i,j](r)*g[i,j](r)*r^2, eps(Float64), 20)
  end

  for i in 1:N, j in 1:N
    if i > j
      I[i,j] = I[j,i]
    end
  end

  return 2π*ρ*sum(I)
end
