"""
perturbation.jl

General functions for perturbation
"""

function pairpotential_minimum(pert::Perturbation, u::Matrix{Pairpotential}, rmin::Matrix{Float64})::Matrix{Float64}
  N::Int = ncomp(pert)

  return [ u[i,j](rmin[i,j]) for i in 1:N, j in 1:N ]
end

function pairpotential_minimum(pert::Perturbation)
  u::Matrix{Pairpotential} = pairpotential(pert)
  rmin::Matrix{Float64} = pairpotential_minimizer(pert)

  return pairpotential_minimum(pert, u, rmin)
end

function hsdiameter_estimate(pert::Perturbation, u::Matrix{Pairpotential}, rmin::Matrix{Float64}, T::Float64)::Matrix{Float64}
  N::Int = ncomp(pert)
  rcore::Matrix{Float64} = coreradius(pert)

  ret = Matrix{Float64}(N,N)
  for i in 1:N, j in 1:N
    i > j && continue
    fopt(σ::Float64)::Float64 = abs(u[i,j](σ) - u[i,j](rmin[i,j]) - kB*T)
    opt = Optim.optimize(fopt, rcore[i,j], rmin[i,j])
    ret[i,j] = ret[j,i] = Optim.minimizer(opt)
  end

  return ret
end

function hsdiameter_estimate(pert::Perturbation, T::Float64)::Matrix{Float64}
  u::Matrix{Pairpotential} = pairpotential(pert)
  rmin::Matrix{Float64} = pairpotential_minimizer(pert)

  return hsdiameter_estimate(pert, u, rmin, T)
end
