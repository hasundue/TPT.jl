"""
nfetb.jl

NFE-TB Perturbation
"""

immutable NFETB{Tn <: NFEPerturbation, Tt <: TBPerturbation} <: Perturbation
  nfe::Tn # NFE-like perturbation
  tb::Tt # TB-like perturbation
end

ncomp(nfetb::NFETB) = ncomp(nfetb.nfe)
coreradius(nfetb::NFETB) = coreradius(nfetb.nfe)
cutoffradius(nfetb::NFETB) = cutoffradius(nfetb.nfe)
hsdiameter_estimate(nfetb::NFETB) = hsdiameter_estimate(nfetb.nfe)

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

function pairpotential_minimizer(nfetb::NFETB)::Array{Float64,2}
  N::Int = ncomp(nfetb)
  rc::Vector{Float64} = nfetb.nfe.pseudo.rc
  u::Array{Function,2} = pairpotential(nfetb)

  ret = Array{Float64,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    r̄c = (rc[i] + rc[j]) / 2

    opt = Optim.optimize(u[i,j], 3*√r̄c, 7*√r̄c)
    ret[i,j] = ret[j,i] =  Optim.minimizer(opt)
  end

  return ret
end

function pairpotential_derivative(nfetb::NFETB)::Array{Function,2}
  N::Int = ncomp(nfetb)
  u′_nfe::Array{Function,2} = pairpotential_derivative(nfetb.nfe)
  u′_tb::Array{Function,2} = pairpotential_derivative(nfetb.tb)

  [ u′(r) = u′_nfe[i,j](r) + u′_tb[i,j](r) for i in 1:N, j in 1:N ]
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
