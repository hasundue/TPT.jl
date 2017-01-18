"""
abstractwca.jl

Common datatypes and functions for WCA and LWCA

References:
"""

abstract AbstractWCA <: DependentReferenceSystem
abstract AbstractOptimizedWCA <: AbstractWCA

ncomp(wca::AbstractWCA)::Int = ncomp(wca.trial)
numberdensity(wca::AbstractWCA)::Vector{Float64} = numberdensity(wca.trial)
totalnumberdensity(wca::AbstractWCA)::Float64 = totalnumberdensity(wca.trial)
composition(wca::AbstractWCA)::Vector{Float64} = composition(wca.trial)
temperature(wca::AbstractWCA)::Float64 = wca.temp
hsdiameter(wca::AbstractWCA)::Array{Float64,2} = hsdiameter(wca.trial)

function emptyradius(wca::AbstractOptimizedWCA)::Array{Float64,2}
  N::Int = ncomp(wca)
  σ::Matrix{Float64} = hsdiameter(wca)
  B::Matrix{Function} = blipfunction(wca)

  ret = Matrix{Float64}(N,N)
  for i in 1:N, j in 1:N
    i > j && continue

    opt = NLopt.Opt(:LN_BOBYQA, 1)
    NLopt.min_objective!(opt, (r, g) -> abs(B[i,j](r[1])))
    NLopt.lower_bounds!(opt, [0.])
    NLopt.upper_bounds!(opt, [σ[i,j]])
    NLopt.stopval!(opt, eps(Float64))

    (fmin, r₀, res) = NLopt.optimize(opt, [σ[i,j]])

    ret[i,j] = ret[j,i] = r₀[1]
  end
  ret
end

# Ref: N. E. Dubinin et al.: Thermochimica Acta, 518 (2011) 9-12.
function entropy(wca::AbstractOptimizedWCA)::Float64
  T::Float64 = temperature(wca)
  S_hs::Float64 = entropy(wca.trial)
  U₀_wca::Float64 = internal(wca)

  ΔS₀_wca = U₀_wca / T
  S = S_hs + ΔS₀_wca
end

# Perturbation-independent part of internal energy
# Ref: N. E. Dubinin et al.: Thermochimica Acta, 518 (2011) 9-12.
function internal(wca::AbstractOptimizedWCA)::Float64
  N::Int = ncomp(wca)
  T::Float64 = temperature(wca)
  ρ::Float64 = totalnumberdensity(wca)
  c::Vector{Float64} = composition(wca)
  hs::IndependentReferenceSystem = wca.trial
  σ::Array{Float64,2} = hsdiameter(hs)
  rmin::Array{Float64,2} = wca.rmin
  rcore::Array{Float64,2} = emptyradius(wca)
  u₀::Array{Function,2} = repulsivepotential(wca)
  g₀_wca::Array{Function,2} = paircorrelation(wca)

  U₀_wca::Float64 = 0

  for i in 1:N, j in 1:N
    i > j && continue

    n = i == j ? 1 : 2
    g₀s_wca::Function = spline(g₀_wca[i,j], rcore[i,j], rmin[i,j], 32)
    u₀s::Function = spline(u₀[i,j], rcore[i,j], rmin[i,j], 32)

    U₀_wca +=
      2π*ρ * n*c[i]*c[j] * ∫(r -> u₀s(r)*g₀s_wca(r)*r^2, rcore[i,j], rmin[i,j])
  end

  return U₀_wca
end
