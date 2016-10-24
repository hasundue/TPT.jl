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

function ncomp(wca::WCASystem)::Int
  ncomp(wca.trial)
end

function ncomp(wca::OptimizedWCASystem)::Int
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

  function fopt(σ::Vector{Float64}, g)::Float64
    ahs = AHSSystem(σ, ρ₀)
    u_hs::Array{Function,2} = pairpotential(ahs)
    y_hs::Array{Function,2} = cavityfunction(ahs)

    for i in 1:N
      B(r) = y_hs[i,i](r) * (exp(-β*u₀t[i,i](r)) - exp(-β*u_hs[i,i](r)))
      I[i] = ∫(B, 0.5σ[i], rmin[i,i])
    end

    return norm(I, 1)
  end

  σ₀d::Vector{Float64} = [σ₀[i,i] for i in 1:N]
  rmind::Vector{Float64} = [rmin[i,i] for i in 1:N]

  opt = Opt(:LN_BOBYQA, N)
  min_objective!(opt, fopt)
  lower_bounds!(opt, 0.5σ₀d)
  upper_bounds!(opt, rmind)
  initial_step!(opt, 0.1rmind)
  xtol_rel!(opt, 1e-3)

  (fmin, σ_wca, ret) = optimize(opt, 0.8rmind)

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

  σ₀::Array{Float64,2} = hsdiameter(wca.trial)
  u₀::Array{Function,2} = wca.u₀
  g_hs::Array{Function,2} = prdf(wca.trial)
  y_hs::Array{Function,2} = cavityfunction(wca.trial)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    u₀t = spline(u₀[i,j], 0.5σ₀[i,j], R_MAX, 256)
    function g(r)
      val = y_hs[i,j](r) * exp(-β*u₀t(r))
      return abs(val) < eps(Float64) ? 0. : val
    end
    ret[i,j] = g
  end

  return ret
end

function psf(wca::OptimizedWCASystem)::Array{Function,2}
  N::Int = ncomp(wca)
  c::Vector{Float64} = composition(wca)

  Sref::Array{Function,2} = psf(wca.trial)
  ρ::Float64 = numberdensity(wca)
  b::Array{Function,2} = blipfunction(wca)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue
    bt = spline(b[i,j], R_MIN, wca.rmin[i,j], 64)
    B(q) = ∫(r -> bt(r) * sin(r*q) / (r*q) * r^2, R_MIN, wca.rmin[i,j], e=1e-3)
    S(q) = Sref[i,j](q) + 4π*ρ * √(c[i]*c[j]) * B(q)
    ret[i,j] = S
  end

  for i in 1:N, j in 1:N
    i < j && continue
    ret[i,j] = ret[j,i]
  end

  return ret
end
