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
  rmin::Array{Float64,2} # position of the minimum of perturbation pair-potential
  u₀::Array{Function,2} # WCA reference pair-potential
  u₁::Array{Function,2} # WCA perturbation pair-potential
end

function ncomp(wca::WCASystem)
  ncomp(wca.trial)
end

function numberdensity(wca::OptimizedWCASystem) :: Float64
  return sum(wca.trial.ρ)
end

function composition(wca::OptimizedWCASystem) :: Vector{Float64}
  ρ = numberdensity(wca)
  return wca.trial.ρ / ρ
end

function TPTSystem(wca::WCASystem{AHSSystem}, pert::Perturbation)
  β = 1 / (kB * wca.T)

  N::Int = ncomp(wca)
  σ₀::Array{Float64,2} = hsdiameter(wca.trial)
  ρ₀::Vector{Float64} = wca.trial.ρ

  u::Array{Function,2} = pairpotential(pert)

  ut = Array{Function}(N,N)
  for i in 1:N, j in 1:N
    ut[i,j] = spline(u[i,j], 0.25σ₀[i,j], 1.5σ₀[i,j], 64)
  end

  rmin = Array{Float64}(N,N)
  u₀ = Array{Function}(N,N)
  u₀t = Array{Function}(N,N)
  u₁ = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    opt = Optim.optimize(ut[i,j], 0.5σ₀[i,j], 1.5σ₀[i,j])
    rmin[i,j] = Optim.minimizer(opt)
    umin::Float64 = Optim.minimum(opt)

    u₀[i,j] = r -> r < rmin[i,j] ? u[i,j](r) - umin : 0.0
    u₀t[i,j] = r -> r < rmin[i,j] ? ut[i,j](r) - umin : 0.0
    u₁[i,j] = r -> r < rmin[i,j] ? umin : u[i,j](r)
  end

  for i in 1:N, j in 1:N
    i < j && continue

    rmin[i,j] = rmin[j,i]
    u₀[i,j] = u₀[j,i]
    u₀t[i,j] = u₀t[j,i]
    u₁[i,j] = u₁[j,i]
  end

  I = Vector{Float64}(N)

  function fopt(σ::Vector{Float64}) :: Float64
    ahs = AHSSystem(σ, ρ₀)
    u_hs::Array{Function,2} = pairpotential(ahs)
    y_hs::Array{Function,2} = cavityfunction(ahs)

    for i in 1:N
      B(r) = y_hs[i,i](r) * (exp(-β*u₀t[i,i](r)) - exp(-β*u_hs[i,i](r)))
      I[i] = ∫(B, 0.5σ[i], rmin[i,i])
    end

    return norm(I)
  end

  if N == 1
    res = Optim.optimize(σ -> fopt([σ]), 0.5σ₀[1], rmin[1], abs_tol=1e-3)
  else
    res = Optim.optimize(fopt, σ₀, f_tol=1e-3, show_trace=true)
  end

  !Optim.converged(res) && error("WCA method couldn't converge")

  # the optimized hard-sphere diameters and the hard-sphere system
  σ_wca = N == 1 ? [Optim.minimizer(res)] : Optim.minimizer(res)

  ahs = AHSSystem(σ_wca, ρ₀)
  optwca = OptimizedWCASystem{AHSSystem}(ahs, wca.T, rmin, u₀, u₁)

  return TPTSystem(optwca, pert)
end

function blipfunction(wca::OptimizedWCASystem)
  β = 1 / (kB * wca.T)

  N = ncomp(wca.trial)

  u₀::Array{Function,2} = wca.u₀
  u_hs::Array{Function,2} = pairpotential(wca.trial)
  y_hs::Array{Function,2} = cavityfunction(wca.trial)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    B(r) = y_hs[i,j](r) * (exp(-β*u₀[i,j](r)) - exp(-β*u_hs[i,j](r)))
    ret[i,j] = B
  end

  return ret
end

function prdf(wca::OptimizedWCASystem{AHSSystem})
  β = 1 / (kB * wca.T)

  N::Int = ncomp(wca.trial)

  u₀::Array{Function,2} = wca.u₀
  g_hs::Array{Function,2} = prdf(wca.trial)
  y_hs::Array{Function,2} = cavityfunction(wca.trial)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    function g(r)
      val = y_hs[i,j](r) * exp(-β*u₀[i,j](r))
      return abs(val) < eps(Float64) ? 0. : val
    end
    ret[i,j] = g
  end

  return ret
end
