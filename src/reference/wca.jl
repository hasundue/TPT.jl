"""
wca.jl

WCA reference system
"""

immutable WCASystem{T <: IndependentReferenceSystem} <: DependentReferenceSystem
  trial::T # trial system with initial conditions
  T::Float64 # Temperature
end

function prdf(wca::WCASystem{AHSSystem}, pert::Perturbation)
  β = 1 / (kB * wca.T)
  σ₀ = wca.trial.σ
  ρ₀ = wca.trial.ρ

  u = pairpotential(pert)

  N = size(u, 1)

  ut = Array{Function}(N)
  for i in 1:N
    ut[i] = tablize(u[i,i], 0.25σ₀[i], 1.5σ₀[i], 100)
  end

  rmin = Vector{Float64}(N)
  u₀ = Vector{Function}(N)
  u₀t = Vector{Function}(N)
  u₁ = Vector{Function}(N)

  for i in 1:N
    opt = Optim.optimize(ut[i], σ₀[i], 1.5σ₀[i])
    rmin[i] = Optim.minimizer(opt)
    umin = Optim.minimum(opt)

    u₀[i] = r -> r < rmin[i] ? u[i,i](r) - umin : 0.0
    u₀t[i] = r -> r < rmin[i] ? ut[i](r) - umin : 0.0
    u₁[i] = r -> r < rmin[i] ? umin : u[i,i](r)
  end

  B = Vector{Function}(N)
  I = Vector{Float64}(N)

  function fopt(σ::Vector{Float64}) :: Float64
    ahs = AHSSystem(σ, ρ₀)
    u_hs = pairpotential(ahs)

    for i in 1:N
      Δf(r) = exp(-β*u₀t[i](r)) - exp(-β*u_hs[i,i](r))
      I[i] = ∫(Δf, 0.5σ[i], rmin[i])
    end

    return norm(I)
  end

  if N == 1
    res = Optim.optimize(σ -> fopt([σ]), 0.5σ₀[1], rmin[1], abs_tol=1e-6)
  else
    res = Optim.optimize(fopt, σ₀, ftol=1e-6)
  end

  if !Optim.converged(res)
    error("WCA method couldn't converge")
  end

  # the optimized hard-sphere diameters and the hard-sphere system
  σ_wca = N == 1 ? [Optim.minimizer(res)] : Optim.minimizer(res)
  ahs = AHSSystem(σ_wca, ρ₀)
  g_hs = prdf(ahs)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    if i == j
      ret[i,j] = r -> r < σ_wca[i] ? 0.0 : g_hs[i,i](r) * exp(-β*u₀[i](r))
    else
      # Not implemented yet
      warn("binary WCA is not implemented yet")
      ret[i,j] = g_hs[i,j]
    end
  end

  return ret
end
