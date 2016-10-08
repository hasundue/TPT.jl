"""
wca.jl

WCA reference system
"""

immutable WCASystem{T <: IndependentReferenceSystem} <: DependentReferenceSystem
  trial::T # trial system with initial conditions
  T::Float64 # Temperature
end

immutable OptimizedWCASystem{T <: IndependentReferenceSystem} <: DependentReferenceSystem
  trial::T # optimized trial system
  T::Float64 # Temperature
  rmin::Vector{Float64} # position of the minimum of perturbation pair-potential
  u₀::Vector{Function} # WCA reference pair-potential
  u₁::Vector{Function} # WCA perturbation pair-potential
end

function TPTSystem(wca::WCASystem{AHSSystem}, pert::Perturbation)
  β = 1 / (kB * wca.T)
  σ₀ = wca.trial.σ
  ρ₀ = wca.trial.ρ

  u = pairpotential(pert)

  N = size(u, 1)

  if N > 1
    warn("WCA with a multi-component AHS system is not implemented yet")
  end

  ut = Array{Function}(N)
  for i in 1:N
    ut[i] = spline(u[i,i], 0.25σ₀[i], 1.5σ₀[i], 64)
  end

  rmin = Vector{Float64}(N)
  u₀ = Vector{Function}(N)
  u₀t = Vector{Function}(N)
  u₁ = Vector{Function}(N)

  for i in 1:N
    opt = Optim.optimize(ut[i], 0.5σ₀[i], 1.5σ₀[i])
    rmin[i] = Optim.minimizer(opt)
    umin = Optim.minimum(opt)

    u₀[i] = r -> r < rmin[i] ? u[i,i](r) - umin : 0.0
    u₀t[i] = r -> r < rmin[i] ? ut[i](r) - umin : 0.0
    u₁[i] = r -> r < rmin[i] ? umin : u[i,i](r)
  end

  I = Vector{Float64}(N)

  function fopt(σ::Vector{Float64}) :: Float64
    ahs = AHSSystem(σ, ρ₀)
    u_hs = pairpotential(ahs)
    y_hs = cavityfunction(ahs)

    for i in 1:N
      y = spline(y_hs[i,i], 0.5σ[i], rmin[i], 64)
      B(r) = y(r) * (exp(-β*u₀t[i](r)) - exp(-β*u_hs[i,i](r)))
      I[i] = ∫(B, 0.5σ[i], rmin[i])
    end

    return norm(I)
  end

  if N == 1
    res = Optim.optimize(σ -> fopt([σ]), 0.5σ₀[1], rmin[1], abs_tol=1e-3)
  else
    res = Optim.optimize(fopt, σ₀, ftol=1e-3)
  end

  if !Optim.converged(res)
    error("WCA method couldn't converge")
  end

  # the optimized hard-sphere diameters and the hard-sphere system
  σ_wca = N == 1 ? [Optim.minimizer(res)] : Optim.minimizer(res)

  ahs = AHSSystem(σ_wca, ρ₀)
  optwca = OptimizedWCASystem{AHSSystem}(ahs, wca.T, rmin, u₀, u₁)

  return TPTSystem(optwca, pert)
end

function blipfunction(wca::OptimizedWCASystem)
  β = 1 / (kB * wca.T)
  u₀ = wca.u₀
  u_hs = pairpotential(wca.trial)
  y_hs = cavityfunction(wca.trial)

  N = length(wca.trial.σ)
  B = Vector{Function}(N)

  for i in 1:N
    B[i] = r -> y_hs[i,i](r) * (exp(-β*u₀[i](r)) - exp(-β*u_hs[i,i](r)))
  end

  return B
end

function prdf(wca::OptimizedWCASystem{AHSSystem})
  β = 1 / (kB * wca.T)
  σ_wca = wca.trial.σ
  u₀ = wca.u₀
  g_hs = prdf(wca.trial)
  y_hs = cavityfunction(wca.trial)

  N = length(σ_wca)
  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    if i == j
      ret[i,j] = r -> y_hs[i,j](r) * exp(-β*u₀[i](r))
    else
      # Not implemented yet
      ret[i,j] = g_hs[i,j]
    end
  end

  return ret
end
