"""
wca.jl

WCA reference system

References:
"""

immutable WCA{T <: IndependentReferenceSystem} <: AbstractWCA
  trial::T # initial guess for trial system
  temp::Float64 # temperature
end

WCA(trial::IndependentReferenceSystem, temp::Real) =
  WCA(trial, convert(Float64, temp))

immutable OptimizedWCA{T <: IndependentReferenceSystem} <: AbstractOptimizedWCA
  trial::T # optimized trial system
  temp::Float64 # temperature
  rmin::Array{Float64,2} # positionof the minimum of pair-potential
  repu::Array{Function,2} # repulsive part of pair-potential
  residue::Float64 # residue in optimization
end

repulsivepotential(wca::OptimizedWCA)::Array{Function,2} = wca.repu

function TPTSystem(wca::WCA{AHS}, pert::Perturbation; kwargs...)
  T::Float64 = temperature(wca)
  β::Float64 = 1 / (kB * T)

  N::Int = ncomp(wca)
  approx::AbstractString = wca.trial.approx
  c::Vector{Float64} = composition(wca)
  σ₀::Matrix{Float64} = hsdiameter(wca.trial)
  ρ₀::Vector{Float64} = numberdensity(wca)

  u::Matrix{Function} = pairpotential(pert)
  rcore::Matrix{Float64} = coreradius(pert)
  rcut::Matrix{Float64} = cutoffradius(pert)
  rmin::Matrix{Float64} = pairpotential_minimizer(pert)
  umin::Matrix{Float64} = pairpotential_minimum(pert, u, rmin)

  # Repulsive-part of pairpotential
  u₀::Array{Function,2} = [ r -> ( r ≤ rmin[i,j] ? u[i,j](r) - umin[i,j] : 0. )
                            for i in 1:N, j in 1:N ]

  I = Vector{Float64}(N)

  function fopt(σ::Vector{Float64}, g)::Float64
    ahs = AHS(approx, σ::Vector, ρ₀::Vector)
    u_hs::Array{Function,2} = pairpotential(ahs)
    g_hs::Array{Function,2} = paircorrelation(ahs)

    for i in 1:N
      y_hs::Function = spline(g_hs[i,i], σ[i], rmin[i,i], 8, bc="extrapolate")
      us::Function = spline(u₀[i,i], rcore[i,i], rmin[i,i], 8)
      B(r) = y_hs(r) * (exp(-β*us(r)) - exp(-β*u_hs[i,i](r)))

      I[i] = ∫(r -> B(r)*r^2, rcore[i,i], rmin[i,i])
    end

    return norm(I, 1)
  end

  opt = NLopt.Opt(:LN_BOBYQA, N)
  NLopt.min_objective!(opt, fopt)
  NLopt.xtol_rel!(opt, 1e-5)
  NLopt.lower_bounds!(opt, [ rcore[i,i] for i in 1:N ])
  NLopt.upper_bounds!(opt, [ 0.99rmin[i,i] for i in 1:N ])
  NLopt.initial_step!(opt, [ 0.01rmin[i,i] for i in 1:N ])

  σ_init = [ hsdiameter_estimate(pert, u, rmin, T)[i,i] for i in 1:N ]
  (fmin, σ_wca, ret) = NLopt.optimize(opt, σ_init)

  ahs = AHS(approx, σ_wca::Vector, ρ₀::Vector)
  optwca = OptimizedWCA(ahs, T, rmin, u₀, fmin)

  TPTSystem(optwca, pert; kwargs...)
end

function blipfunction(wca::OptimizedWCA)::Array{Function,2}
  N = ncomp(wca.trial)

  T = temperature(wca)
  β = 1 / (kB * T)

  σ::Array{Float64,2} = hsdiameter(wca.trial)
  rmin::Array{Float64,2} = wca.rmin

  u₀::Array{Function,2} = repulsivepotential(wca)
  u_hs::Array{Function,2} = pairpotential(wca.trial)
  g_hs::Array{Function,2} = paircorrelation(wca.trial)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue
    y_hs::Function = spline(g_hs[i,j], σ[i,j], rmin[i,j], 8, bc="extrapolate")
    B(r) = y_hs(r) * (exp(-β*u₀[i,j](r)) - exp(-β*u_hs[i,j](r)))
    ret[i,j] = ret[j,i] = B
  end

  return ret
end

function paircorrelation(wca::OptimizedWCA)
  T = temperature(wca)
  β = 1 / (kB * T)

  N::Int = ncomp(wca.trial)

  σ₀::Array{Float64,2} = hsdiameter(wca)
  rmin::Array{Float64,2} = wca.rmin
  u₀::Array{Function,2} = repulsivepotential(wca)
  g_hs::Array{Function,2} = paircorrelation(wca.trial)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    rcut = 0.5σ₀[i,j]

    y_hs = spline(g_hs[i,j], σ₀[i,j], rmin[i,j], 8, bc="extrapolate")

    function g(r)::Float64
      if r < rcut
        0
      elseif r < σ₀[i,j]
        y_hs(r) * exp(-β*u₀[i,j](r))
      elseif r < rmin[i,j]
        g_hs[i,j](r) * exp(-β*u₀[i,j](r))
      else
        g_hs[i,j](r)
      end
    end

    ret[i,j] = ret[j,i] = g
  end

  return ret
end

function structurefactor(wca::OptimizedWCA)::Array{Function,2}
  N::Int = ncomp(wca)
  c::Vector{Float64} = composition(wca)

  S₀::Array{Function,2} = structurefactor(wca.trial)
  σ::Array{Float64,2} = hsdiameter(wca.trial)
  ρ::Float64 = totalnumberdensity(wca)
  rcore::Array{Float64,2} = emptyradius(wca)
  rmin::Array{Float64,2} = wca.rmin
  B::Array{Function,2} = blipfunction(wca)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    ϵ = eps(Float64)

    Bs1 = spline(B[i,j], rcore[i,j], σ[i,j]-ϵ, 8)
    Bs2 = spline(B[i,j], σ[i,j]+ϵ, rmin[i,j], 8)

    B̃(q) =
      4π * ∫(r -> Bs1(r) * sin(r*q) / (r*q) * r^2, rcore[i,j], σ[i,j]-ϵ) +
      4π * ∫(r -> Bs2(r) * sin(r*q) / (r*q) * r^2, σ[i,j]+ϵ, rmin[i,j])

    # S(q) = S₀[i,j](q) + ρ * √(c[i]*c[j]) * B̃(q)
    S(q) = S₀[i,j](q) / (1 - ρ * √(c[i]*c[j]) * S₀[i,j](q) * B̃(q))

    ret[i,j] = ret[j,i] = S
  end

  return ret
end
