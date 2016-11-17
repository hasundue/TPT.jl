"""
wca.jl

WCA reference system
"""

abstract AbstractWCA <: DependentReferenceSystem

immutable WCA{T <: IndependentReferenceSystem} <: AbstractWCA
  trial::T # initial guess for trial system
  temp::Float64 # temperature
end

WCA(trial::IndependentReferenceSystem, temp::Real) = WCA(trial, convert(Float64, temp))

immutable OptimizedWCA{T <: IndependentReferenceSystem} <: AbstractWCA
  trial::T # optimized trial system
  temp::Float64 # temperature
  rmin::Array{Float64,2} # positionof the minimum of pair-potential
  repu::Array{Function,2} # repulsive part of pair-potential
  tail::Array{Function,2} # tail part of pair-potential
end

ncomp(wca::AbstractWCA)::Int = ncomp(wca.trial)
numberdensity(wca::AbstractWCA)::Vector{Float64} = numberdensity(wca.trial)
totalnumberdensity(wca::AbstractWCA)::Float64 = totalnumberdensity(wca.trial)
composition(wca::AbstractWCA)::Vector{Float64} = composition(wca.trial)
temperature(wca::AbstractWCA)::Float64 = wca.temp
hsdiameter(wca::AbstractWCA)::Array{Float64,2} = hsdiameter(wca.trial)

repulsivepotential(wca::OptimizedWCA)::Array{Function,2} = wca.repu
tailpotential(wca::AbstractWCA)::Array{Function,2} = wca.tail

function TPTSystem(wca::WCA{AHS}, pert::Perturbation; kwargs...)
  T = temperature(wca)
  β = 1 / (kB * T)

  N::Int = ncomp(wca)
  σ₀::Array{Float64,2} = hsdiameter(wca.trial)
  ρ₀::Vector{Float64} = numberdensity(wca)

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
    rmin[i,j] = rmin[j,i] =  Optim.minimizer(opt)
    umin::Float64 = Optim.minimum(opt)

    u₀[i,j] = u₀[j,i] = r -> r < rmin[i,j] ? u[i,j](r) - umin : 0.0
    u₀t[i,j] = u₀t[j,i] = r -> r < rmin[i,j] ? ut[i,j](r) - umin : 0.0
    u₁[i,j] = u₁[i,j] = r -> r < rmin[i,j] ? umin : u[i,j](r)
  end

  I = Vector{Float64}(N)

  function fopt(σ::Vector{Float64}, g)::Float64
    ahs = AHS(σ, ρ₀)
    u_hs::Array{Function,2} = pairpotential(ahs)
    y_hs::Array{Function,2} = cavityfunction(ahs)

    for i in 1:N
      B(r) = y_hs[i,i](r) * (exp(-β*u₀t[i,i](r)) - exp(-β*u_hs[i,i](r)))
      I[i] = ∫(r -> B(r)*r^2, 0.5σ[i], rmin[i,i])
    end

    return norm(I, 1)
  end

  σ₀d::Vector{Float64} = [σ₀[i,i] for i in 1:N]
  rmind::Vector{Float64} = [rmin[i,i] for i in 1:N]
  σ_init::Vector{Float64} = [min(σ₀d[i], rmind[i]) for i in 1:N]

  opt = NLopt.Opt(:LN_BOBYQA, N)
  NLopt.min_objective!(opt, fopt)
  NLopt.lower_bounds!(opt, 0.8σ_init)
  NLopt.upper_bounds!(opt, rmind)
  NLopt.initial_step!(opt, 0.01σ_init)
  NLopt.xtol_rel!(opt, 1e-3)

  (fmin, σ_wca, ret) = NLopt.optimize(opt, σ_init)

  ahs = AHS(σ_wca, ρ₀)
  optwca = OptimizedWCA(ahs, T, rmin, u₀, u₁)

  TPTSystem(optwca, pert; kwargs...)
end

function blipfunction(wca::OptimizedWCA)
  T = temperature(wca)
  β = 1 / (kB * T)

  N = ncomp(wca.trial)

  u₀::Array{Function,2} = repulsivepotential(wca)
  u_hs::Array{Function,2} = pairpotential(wca.trial)
  y_hs::Array{Function,2} = cavityfunction(wca.trial)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    B(r) = y_hs[i,j](r) * (exp(-β*u₀[i,j](r)) - exp(-β*u_hs[i,j](r)))

    ret[i,j] = ret[j,i] =  B
  end

  return ret
end

function paircorrelation(wca::OptimizedWCA)
  β = 1 / (kB * wca.temp)

  N::Int = ncomp(wca.trial)

  σ₀::Array{Float64,2} = hsdiameter(wca.trial)
  u₀::Array{Function,2} = repulsivepotential(wca)
  g_hs::Array{Function,2} = paircorrelation(wca.trial)
  y_hs::Array{Function,2} = cavityfunction(wca.trial)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    u₀t = spline(u₀[i,j], 0.5σ₀[i,j], R_MAX, 256)

    function g(r)
      val = y_hs[i,j](r) * exp(-β*u₀t(r))
      abs(val) < eps(Float64) ? 0. : val
    end

    ret[i,j] = ret[j,i] = g
  end

  return ret
end

function structurefactor(wca::WCA)::Array{Function,2}
  N::Int = ncomp(wca)
  c::Vector{Float64} = composition(wca)

  Sref::Array{Function,2} = structurefactor(wca.trial)
  ρ::Float64 = numberdensity(wca)
  b::Array{Function,2} = blipfunction(wca)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    bt = spline(b[i,j], R_MIN, wca.rmin[i,j], 64)

    B(q) = ∫(r -> bt(r) * sin(r*q) / (r*q) * r^2, R_MIN, wca.rmin[i,j], e=1e-3)
    S(q) = Sref[i,j](q) / (1 - 4π*ρ * √(c[i]*c[j]) * Sref[i,j](q) * B(q))

    ret[i,j] = ret[j,i] = S
  end

  return ret
end

# Ref: N. E. Dubinin et al.: Thermochimica Acta, 518 (2011) 9-12.
function entropy(wca::OptimizedWCA)::Float64
  T::Float64 = temperature(wca)
  S_hs::Float64 = entropy(wca.trial)
  U₀_wca::Float64 = internal(wca)

  ΔS₀_wca = U₀_wca / T

  S = S_hs + ΔS₀_wca
end

# Perturbation-independent part of internal energy
# Ref: N. E. Dubinin et al.: Thermochimica Acta, 518 (2011) 9-12.
function internal(wca::OptimizedWCA)::Float64
  N::Int = ncomp(wca)
  T::Float64 = temperature(wca)
  ρ::Float64 = totalnumberdensity(wca)
  c::Vector{Float64} = composition(wca)
  hs::IndependentReferenceSystem = wca.trial
  σ::Array{Float64,2} = hsdiameter(hs)
  rmin::Array{Float64,2} = wca.rmin
  u₀::Array{Function,2} = repulsivepotential(wca)
  g₀_wca::Array{Function,2} = paircorrelation(wca)

  U₀_wca::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue

    a = i == j ? 1 : 2
    u₀t = spline(u₀[i,j], σ[i,j]/2, rmin[i,j], 32)

    U₀_wca += a * 2π*ρ * c[i]*c[j] * ∫(r -> u₀t(r)*g₀_wca[i,j](r)*r^2, σ[i,j]/2, rmin[i,j])
  end

  U₀_wca
end

# Perturbation-dependent part of internal energy
# Ref: N. E. Dubinin et al.: Thermochimica Acta, 518 (2011) 9-12.
function internal(wca::OptimizedWCA, pert::Perturbation)::Float64
  N::Int = ncomp(wca)
  c::Vector{Float64} = composition(wca)
  hs::IndependentReferenceSystem = wca.trial
  σ::Array{Float64,2} = hsdiameter(hs)
  rmin::Array{Float64,2} = wca.rmin
  ρ::Float64 = totalnumberdensity(wca)
  u::Array{Function,2} = pairpotential(pert)
  u₀::Array{Function,2} = repulsivepotential(wca)
  g_HS::Array{Function,2} = paircorrelation(hs)

  U::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue

    a = i == j ? 1 : 2
    A = a * 2π*ρ * c[i]*c[j]

    ut = spline(u[i,j], σ[i,j], R_MAX, 256)
    u₀t = spline(u₀[i,j], σ[i,j], rmin[i,j], 16)

    U += A * ∫(r -> ut(r)*g_HS[i,j](r)*r^2, σ[i,j], R_MAX)
    U -= A * ∫(r -> u₀t(r)*g_HS[i,j](r)*r^2, σ[i,j], rmin[i,j])
  end

  return U
end
