#
# Common functions for ReferenceSystem
#
function nndistance(ref::ReferenceSystem, g::Array{Function,2})::Array{Float64,2}
  N::Int = ncomp(ref)
  σ::Array{Float64,2} = hsdiameter(ref)
  d = Array{Float64,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    opt = Optim.optimize(r -> -g[i,j](r), 0.5σ[i,j], 1.5σ[i,j])
    d[i,j] = d[j,i] = Optim.minimizer(opt)
  end

  return d
end

function nndistance(ref::ReferenceSystem)::Array{Float64,2}
  g::Array{Function,2} = paircorrelation(ref)
  return nndistance(ref, g)
end

# pairwise contribution to internal energy
function internal(ref::ReferenceSystem, pert::Perturbation)::Float64
  N::Int = ncomp(ref)
  ρ::Float64 = totalnumberdensity(ref)
  c::Vector{Float64} = composition(ref)
  g::Array{Function,2} = paircorrelation(ref)
  rcore::Array{Float64,2} = emptyradius(ref)

  u::Array{Function,2} = pairpotential(pert)
  rcut::Array{Float64,2} = cutoffradius(pert)

  U::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue
    n = i == j ? 1 : 2

    gs = spline(g[i,j], rcore[i,j], rcut[i,j], 32)
    us = spline(u[i,j], rcore[i,j], rcut[i,j], 32)
    U += 2π*ρ * n*c[i]*c[j] * ∫(r -> us(r)*gs(r)*r^2, rcore[i,j], rcut[i,j])
  end

  return U
end
