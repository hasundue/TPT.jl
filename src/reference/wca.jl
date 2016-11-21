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
  T::Float64 = temperature(wca)
  β::Float64 = 1 / (kB * T)

  N::Int = ncomp(wca)
  σ₀::Array{Float64,2} = hsdiameter(wca.trial)
  ρ₀::Vector{Float64} = numberdensity(wca)

  u::Array{Function,2} = pairpotential(pert)

  us = Array{Function,2}(N,N)
  for i in 1:N, j in 1:N
    i > j && continue
    us[i,j] = us[j,i] = spline(u[i,j], 0.25σ₀[i,j], 1.5σ₀[i,j], 32)
  end

  rmin = Array{Float64,2}(N,N)
  u₀ = Array{Function,2}(N,N)
  u₀s = Array{Function,2}(N,N)
  u₁ = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    opt = Optim.optimize(us[i,j], 0.5σ₀[i,j], 1.5σ₀[i,j])
    rmin[i,j] = rmin[j,i] =  Optim.minimizer(opt)
    umin::Float64 = Optim.minimum(opt)

    u₀[i,j] = u₀[j,i] = r -> r < rmin[i,j] ? u[i,j](r) - umin : 0.0
    u₀s[i,j] = u₀s[j,i] = r -> r < rmin[i,j] ? us[i,j](r) - umin : 0.0
    u₁[i,j] = u₁[i,j] = r -> r < rmin[i,j] ? umin : u[i,j](r)
  end

  I = Vector{Float64}(N)

  function fopt(σ::Vector{Float64}, g)::Float64
    ahs = AHS(σ::Vector, ρ₀::Vector)
    u_hs::Array{Function,2} = pairpotential(ahs)
    g_hs::Array{Function,2} = paircorrelation(ahs)

    ϵ = eps(Float64)

    for i in 1:N
      rcut::Float64 = 0.5σ[i]
      y_hs::Function = spline(g_hs[i,i], σ[i], rmin[i,i], 8, bc="extrapolate")
      B(r) = y_hs(r) * (exp(-β*u₀s[i,i](r)) - exp(-β*u_hs[i,i](r)))

      I[i] = ∫(r -> B(r)*r^2, rcut, σ[i]-ϵ) +
             ∫(r -> B(r)*r^2, σ[i]+ϵ, rmin[i,i])
    end

    return norm(I, 1)
  end

  σ₀d::Vector{Float64} = [ σ₀[i,i] for i in 1:N ]
  rmind::Vector{Float64} = [ rmin[i,i] for i in 1:N ]
  σ_init::Vector{Float64} = [ min(σ₀d[i], rmind[i]) for i in 1:N ]

  if N == 1
    opt = NLopt.Opt(:GN_DIRECT, N)
  else
    opt = NLopt.Opt(:LN_BOBYQA, N)
  end

  NLopt.min_objective!(opt, fopt)
  NLopt.lower_bounds!(opt, 0.5σ₀d)
  NLopt.upper_bounds!(opt, rmind)
  NLopt.stopval!(opt, 1e-1)

  (fmin, σ_wca, ret) = NLopt.optimize(opt, σ_init)

  ahs = AHS(σ_wca::Vector, ρ₀::Vector)
  optwca = OptimizedWCA(ahs, T, rmin, u₀, u₁)

  TPTSystem(optwca, pert; kwargs...)
end

function blipfunction(wca::OptimizedWCA)
  T = temperature(wca)
  β = 1 / (kB * T)

  N = ncomp(wca.trial)
  σ::Array{Float64,2} = hsdiameter(wca.trial)
  rmin::Array{Float64,2} = wca.rmin

  u₀::Array{Function,2} = repulsivepotential(wca)
  u_hs::Array{Function,2} = pairpotential(wca.trial)
  g_hs::Array{Function,2} = paircorrelation(wca.trial)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue
    y_hs::Function = spline(g_hs[i,j], σ[i,i], rmin[i,i], 8, bc="extrapolate")
    B(r) = y_hs(r) * (exp(-β*u₀[i,j](r)) - exp(-β*u_hs[i,j](r)))
    ret[i,j] = ret[j,i] =  B
  end

  return ret
end

function paircorrelation(wca::OptimizedWCA)
  T = temperature(wca)
  β = 1 / (kB * T)

  N::Int = ncomp(wca.trial)

  σ₀::Array{Float64,2} = hsdiameter(wca.trial)
  rmin::Array{Float64,2} = wca.rmin
  u₀::Array{Function,2} = repulsivepotential(wca)
  g_hs::Array{Function,2} = paircorrelation(wca.trial)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    rcut = 0.5σ₀[i,j]

    y_hs = spline(g_hs[i,j], σ₀[i,j], rmin[i,j], 4, bc="extrapolate")
    u₀s = spline(u₀[i,j], rcut, rmin[i,j], 16)

    function g(r)::Float64
      if r < rcut
        0
      elseif r < rmin[i,j]
        y_hs(r) * exp(-β*u₀s(r))
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

  Sref::Array{Function,2} = structurefactor(wca.trial)
  σ::Array{Float64,2} = hsdiameter(wca.trial)
  ρ::Float64 = totalnumberdensity(wca)
  rmin::Array{Float64,2} = wca.rmin
  b::Array{Function,2} = blipfunction(wca)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    ϵ = eps(Float64)
    rcut = σ[i,j]/2

    bs1 = spline(b[i,j], rcut, σ[i,j]-ϵ, 4)
    bs2 = spline(b[i,j], σ[i,j]+ϵ, rmin[i,j], 4)

    B(q) =
      ∫(r -> bs1(r) * sin(r*q) / (r*q) * r^2, rcut, σ[i,j]-ϵ, e=1e-3) +
      ∫(r -> bs2(r) * sin(r*q) / (r*q) * r^2, σ[i,j]+ϵ, rmin[i,j], e=1e-3)

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
    g₀s_wca::Function = spline(g₀_wca[i,j], σ[i,j]/2, rmin[i,j], 16)
    u₀s::Function = spline(u₀[i,j], σ[i,j]/2, rmin[i,j], 16)

    U₀_wca += a * 2π*ρ * c[i]*c[j] * ∫(r -> u₀s(r)*g₀s_wca(r)*r^2, σ[i,j]/2, rmin[i,j])
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
  g_hs::Array{Function,2} = paircorrelation(hs)

  U::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue

    a = i == j ? 1 : 2
    A = a * 2π*ρ * c[i]*c[j]

    us = spline(u[i,j], σ[i,j], R_MAX, 32)
    u₀s = spline(u₀[i,j], σ[i,j], rmin[i,j], 16)
    gs_hs = spline(g_hs[i,j], σ[i,j], R_MAX, 32)

    U += A * ∫(r -> us(r)*gs_hs(r)*r^2, σ[i,j], R_MAX)
    U -= A * ∫(r -> u₀s(r)*gs_hs(r)*r^2, σ[i,j], rmin[i,j])
  end

  return U
end
