"""
lwca.jl

Linearized WCA (LWCA) reference system

References:

* A. Meyer et al. CHemical Physics 49 (1980) 147-152
"""

immutable LWCA{T <: IndependentReferenceSystem} <: AbstractWCA
  trial::T
  temp::Float64
  struct::Symbol
end

LWCA(trial::IndependentReferenceSystem, temp::Real; struct = :linear) =
  LWCA(trial, convert(Float64, temp), struct)

immutable OptimizedLWCA{T <: IndependentReferenceSystem} <: AbstractOptimizedWCA
  trial::T
  temp::Float64
  rmin::Matrix{Float64}
  C₁::Matrix{Float64}
  C₂::Matrix{Float64}
  K₁::Matrix{Float64}
  K₂::Matrix{Float64}
  residue::Float64 # residue in optimization
end

function TPTSystem(lwca::LWCA{AHS}, pert::Perturbation; kwargs...)
  N::Int = ncomp(lwca)

  T::Float64 = temperature(lwca)
  β::Float64 = 1 / (kB * T)

  approx::Symbol = lwca.trial.approx
  c::Vector{Float64} = composition(lwca)
  ρ₀::Vector{Float64} = numberdensity(lwca.trial)

  u::Array{Function,2} = pairpotential(pert)
  u′::Array{Function,2} = pairpotential_derivative(pert)

  rcore::Array{Float64,2} = coreradius(pert)
  rcut::Array{Float64,2} = cutoffradius(pert)
  rmin::Array{Float64,2} = pairpotential_minimizer(pert)
  umin::Array{Float64,2} = [ u[i,j](rmin[i,j]) for i in 1:N, j in 1:N ]

  # Repulsive-part of pairpotential
  u₀::Array{Function,2} = [ r -> ( r < rmin[i,j] ? u[i,j](r) - umin[i,j] : 0. )
                            for i in 1:N, j in 1:N ]

  function fopt(σ::Vector{Float64})::Float64
    ahs = AHS(approx, σ, ρ₀)
    g₀ = contactvalue(ahs)
    g′₀ = contactgradient(ahs)

    residue::Float64 = 0

    for i in 1:N
      Y = g′₀[i,i] / g₀[i,i] * σ[i]

      X = (-2β * σ[i] * u′[i,i](σ[i]) + Y + 2) /
          (-β * σ[i] * u′[i,i](σ[i]) + Y + 2)

      if X ≤ 0
        residue += Inf
      else
        residue += abs(β * u₀[i,i](σ[i]) - log(X))
      end
    end

    return residue
  end

  opt = NLopt.Opt(:GN_DIRECT, N)
  NLopt.min_objective!(opt, (σ, g) -> fopt(σ))
  NLopt.stopval!(opt, 1e-5)
  NLopt.lower_bounds!(opt, [ rcore[i,i] for i in 1:N ])
  NLopt.upper_bounds!(opt, [ 0.99rmin[i,i] for i in 1:N ])
  NLopt.initial_step!(opt, [ 0.05rmin[i,i] for i in 1:N ])

  σ_init = [ hsdiameter_estimate(pert, u, rmin, T)[i,i] for i in 1:N ]
  (fmin, σ_wca, res) = NLopt.optimize(opt, σ_init)

  ahs = AHS(approx, σ_wca, ρ₀)

  σ₀ = hsdiameter(ahs)
  g₀ = contactvalue(ahs)
  g′₀ = contactgradient(ahs)
  Y = [ g′₀[i,j] / g₀[i,j] * σ₀[i,j] for i in 1:N, j in 1:N ]

  C₁ = Array{Float64,2}(N,N) # C-
  C₂ = Array{Float64,2}(N,N) # C+
  K₁ = Array{Float64,2}(N,N) # K-
  K₂ = Array{Float64,2}(N,N) # K+

  for i in 1:N, j in 1:N
    i > j && continue

    σᵢⱼ::Float64 = (σ_wca[i] + σ_wca[j]) / 2
    u₀_σ::Float64 = u₀[i,j](σᵢⱼ)
    u′_σ::Float64 = u′[i,j](σᵢⱼ)

    C₁[i,j] = C₁[j,i] = σᵢⱼ^2 * g₀[i,j] * exp(-β*u₀_σ)
    C₂[i,j] = C₂[j,i] = σᵢⱼ^2 * g₀[i,j] * (exp(-β*u₀_σ) - 1)

    K₁[i,j] = K₁[j,i] = -β*σᵢⱼ*u′_σ + Y[i,j] + 2
    K₂[i,j] = K₂[j,i] =
      -β * σᵢⱼ * u′_σ * exp(-β*u₀_σ) / (exp(-β*u₀_σ) - 1) + Y[i,j] + 2
  end

  if lwca.struct == :linear
    optwca = OptimizedLWCA(ahs, T, rmin, C₁, C₂, K₁, K₂, fmin)
  elseif lwca.struct == :full
    optwca = OptimizedWCA(ahs, T, rmin, u₀, fmin)
  end

  TPTSystem(optwca, pert; kwargs...)
end

function blipfunction(lwca::OptimizedLWCA)::Array{Function,2}
  @attach(lwca, C₁, C₂, K₁, K₂)

  N = ncomp(lwca.trial)

  T = temperature(lwca)
  β = 1 / (kB * T)

  σ::Array{Float64,2} = hsdiameter(lwca)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    C₁′ = K₁[i,j] * C₁[i,j] / σ[i,j]
    C₂′ = K₂[i,j] * C₂[i,j] / σ[i,j]

    R₁ = σ[i,j] / K₁[i,j]
    R₂ = - σ[i,j] / K₂[i,j]

    function C(r)::Float64
      if r < σ[i,j] - R₁
        0
      elseif r < σ[i,j]
        C₁[i,j] + C₁′ * (r - σ[i,j])
      elseif r < σ[i,j] + R₂
        C₂[i,j] + C₂′ * (r - σ[i,j])
      else
        0
      end
    end

    B(r)::Float64 = C(r) / r^2

    ret[i,j] = ret[j,i] = B
  end

  return ret
end

function paircorrelation(lwca::OptimizedLWCA)::Array{Function,2}
  N::Int = ncomp(lwca)

  g_hs::Array{Function,2} = paircorrelation(lwca.trial)
  B::Array{Function,2} = blipfunction(lwca)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    g(r)::Float64 = g_hs[i,j](r) + B[i,j](r)

    ret[i,j] = ret[j,i] = g
  end

  return ret
end

function structurefactor(lwca::OptimizedLWCA)::Array{Function,2}
  @attach(lwca, C₁, C₂, K₁, K₂)

  N::Int = ncomp(lwca)
  c::Vector{Float64} = composition(lwca)

  S₀::Array{Function,2} = structurefactor(lwca.trial)
  ρ::Float64 = totalnumberdensity(lwca)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    σ::Float64 = hsdiameter(lwca)[i,j]

    B̃(q)::Float64 = -4π / factorial(3) * q*σ * sin(q*σ) / (q*σ) *
                    ( σ*C₁[i,j] / K₁[i,j]^2 - σ*C₂[i,j] / K₂[i,j]^2 )

    # S(q) = S₀[i,j](q) + ρ * √(c[i]*c[j]) * B̃(q)
    S(q) = S₀[i,j](q) / (1 - ρ * √(c[i]*c[j]) * S₀[i,j](q) * B̃(q))

    ret[i,j] = ret[j,i] = S
  end

  return ret
end
