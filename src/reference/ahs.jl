"""
ahs.jl

Types and functions for additive hard-sphere reference systems.

References:

* S. B. Yuste et al: J. Chem. Phys., Vol. 108, No.9, 1 (1998), 3683-3693.
* K. Hoshino. and W. H. Young: J. Phys. F: Metal Phys., Vol. 10 (1980), 1365-1374.

"""

immutable AHSSystem <: IndependentReferenceSystem
  σ::Vector{Float64} # hard-sphere diameters
  ρ::Vector{Float64} # number densities
end # type

@doc doc"""
AHSSystem(; kwargs...)

Initializes an additive hard-sphere reference system.
You may use one of possible combinations of keyword arguments.

One-component case:

* two out of σ, ρ, and η (all scalar)

Multi-component case:

* σ (vector) and ρ (vector)
* ρ (scalar), σ (vector), and c (vector)
* η (scalar), σ (vector), and c (vector)

You may use sigma, rho, or eta instead of σ, ρ, or η, respectively.
""" ->
function AHSSystem(; kwargs...)
  keys, vals = expandkwargs(kwargs)
  @assert 2 ≤ length(keys) ≤ 3 "wrong number of arguments"

  if all(x -> !(typeof(x) <: Vector), vals) # one-component case

    if greekin(:σ, keys) && greekin(:ρ, keys)
      σ = vals[greekfind(keys, :σ)]
      ρ = vals[greekfind(keys, :ρ)]

    elseif greekin(:ρ, keys) && greekin(:η, keys)
      ρ = vals[greekfind(keys, :ρ)]
      η = vals[greekfind(keys, :η)]
      σ = (6/π * η / ρ)^(1/3)

    elseif greekin(:η, keys) && greekin(:σ, keys)
      η = vals[greekfind(keys, :η)]
      σ = vals[greekfind(keys, :σ)]
      ρ = 6/π * η / σ^3

    else
      error("invalid arguments")
    end # if

    return AHSSystem([σ], [ρ])

  else # Multi-component case

    if length(kwargs) == 2
      @assert all(x-> typeof(x) <: Vector, vals) "invalid arguments"
      σ = vals[greekfind(keys, :σ)]
      ρ = vals[greekfind(keys, :ρ)]

    else
      @assert greekin(:σ, keys) && greekin(:c, keys) "invalid arguments"
      σ = vals[greekfind(keys, :σ)]
      c = vals[findfirst(keys, :c)]

      if greekin(:ρ, keys)
        ρ = c .* vals[greekfind(keys, :ρ)]

      elseif greekin(:η, keys)
        η = vals[greekfind(keys, :η)]
        α = σ / maximum(σ)
        ρ = (6/π) ./ σ.^3 .* (c .* α.^3) ./ sum(c .* α.^3) * η
      end
    end

    @assert length(σ) == length(ρ) "invalid arguments"
    return AHSSystem(σ, ρ)
  end
end

function ncomp(ahs::AHSSystem)::Int
  length(ahs.σ)
end

function composition(ahs::AHSSystem)::Vector{Float64}
  ahs.ρ / sum(ahs.ρ)
end

function numberdensity(ahs::AHSSystem)::Float64
  sum(ahs.ρ)
end

function totalnumberdensity(ahs::AHSSystem)::Float64
  sum(ahs.ρ)
end

function hsdiameter(sys::AHSSystem)::Array{Float64,2}
  N = ncomp(sys)

  ret = Array{Float64}(N,N)

  for i in 1:N, j in 1:N
    if i == j
      ret[i,j] = sys.σ[i]
    else
      ret[i,j] = (sys.σ[i] + sys.σ[j]) / 2
    end
  end

  return ret
end

function packingfraction(ahs::AHSSystem)::Vector{Float64}
  @attach(ahs, ρ, σ)
  return (π/6) * ρ .* σ.^3
end

function totalpackingfraction(ahs::AHSSystem)::Float64
  η = packingfraction(ahs)
  return sum(η)
end

function temperature(ahs::AHSSystem)::Float64
  return ahs.T
end

function pairpotential(sys::AHSSystem)
  σ = sys.σ
  N = length(σ)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    function u(r) :: Float64
      if r < (σ[i] + σ[j]) / 2
        Inf
      else
        0
      end
    end
    ret[i,j] = u
  end

  return ret
end

# ref: S. B. Yuste et al: J. Chem. Phys., Vol. 108, No.9, 1 (1998), 3683-3693.
function lrdf(sys::AHSSystem)
  if length(sys.σ) == 1
    lrdf1(sys)
  elseif length(sys.σ) == 2
    lrdf2(sys)
  else
    error("unsupported number of components")
  end
end

function lrdf1(sys::AHSSystem)
  σ = sys.σ[1]
  ρ = sys.ρ[1]
  η = π/6 * ρ * σ^3

  L⁰ = 2π * (1 + 2η) / (1 - η)^2
  L¹ = 2π*σ * (1 + 1/2*η) / (1 - η)^2

  ϕ₁(x) = x^-2 * (1 - x - exp(-x))
  ϕ₂(x) = x^-3 * (1 - x + x^2/2 - exp(-x))

  G(s) = exp(-σ*s) / (2π*s^2) * (L⁰ + L¹*s) /
  (1 - ρ*(ϕ₂(σ*s)*σ^3*L⁰ + ϕ₁(σ*s)*σ^2*L¹))

  ret = Array{Function,2}(1,1)
  ret[1,1] = G

  return ret
end

function lrdf2(sys::AHSSystem)
  σ = sys.σ
  ρ = sys.ρ

  # auxiliary functions
  L = Array{Function}(2,2)
  A = Array{Function}(2,2)

  # scalar constans
  ζ₁ = sum(ρ .* σ)
  ζ₂ = sum(ρ .* σ.^2)
  ζ₃ = sum(ρ .* σ.^3)
  η = (π/6) * ζ₃
  λ = 2π / (1-η)
  λ′ = π^2 * ζ₂ / (1-η)^2

  # auxiliary functions
  ϕ₁(x) = x^-2 * (1 - x - exp(-x))
  ϕ₂(x) = x^-3 * (1 - x + x^2/2 - exp(-x))

  for i in 1:2, j in 1:2
    σᵢⱼ = (σ[i] + σ[j]) / 2
    L⁰ = λ + λ′*σ[j]
    L¹ = λ*σᵢⱼ + (1/2)*λ′*σ[i]*σ[j]
    L[i,j] = s -> L⁰ + L¹*s
    A[i,j] = s -> ρ[i]*(ϕ₂(σ[i]*s)*σ[i]^3*L⁰ + ϕ₁(σ[i]*s)*σ[i]^2*L¹)
  end

  D(s) = (1 - A[1,1](s)) * (1 - A[2,2](s)) - A[1,2](s) * A[2,1](s)

  σ₁₂ = (σ[1] + σ[2]) / 2

  G₁₁(s) = exp(-σ[1]*s) / (2π*s^2) *
           (L[1,1](s)*(1 - A[2,2](s)) + L[1,2](s)*A[2,1](s)) / D(s)

  G₁₂(s) = exp(-σ₁₂*s) / (2π*s^2) *
           (L[1,2](s)*(1 - A[1,1](s)) + L[1,1](s)*A[1,2](s)) / D(s)

  G₂₁(s) = exp(-σ₁₂*s) / (2π*s^2) *
           (L[2,1](s)*(1 - A[2,2](s)) + L[2,2](s)*A[2,1](s)) / D(s)

  G₂₂(s) = exp(-σ[2]*s) / (2π*s^2) *
           (L[2,2](s)*(1 - A[1,1](s)) + L[2,1](s)*A[1,2](s)) / D(s)

  return [G₁₁ G₁₂; G₂₁ G₂₂]

end


@doc doc"""
psf(sys::AHSSystem)

Returns a structure factor of additive hard-sphere system in wave-number space.
Currently the cases where N = 1 or N = 2 are only supported.
""" ->
function psf(sys::AHSSystem)
  if length(sys.σ) == 1
    return psf1(sys)
  end

  ρ = sys.ρ

  N = length(ρ) # number of components

  # RDF in Laplace space
  if N == 2
    G = lrdf2(sys)
  else
    error("unsupported number of components")
  end

  x = ρ / sum(ρ) # concentrations

  I = eye(N) # identity matrix

  # matrix of structure factors
  F = Array{Function}(N,N)
  h = Array{Function}(N,N)
  S = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    F[i,j] = s -> (s^2 * G[i,j](s) - 1) / s^3
    h[i,j] = q -> -4π * real(F[i,j](im*q))
    # S[i,j] = q -> x[i]*I[i,j] + sum(ρ)*x[i]*x[j]*h[i,j](q) # Faber-Ziman?
    S[i,j] = q -> I[i,j] + √(ρ[i]*ρ[j])*h[i,j](q) # Ashcroft-Langreth
  end

  return S
end

function psf1(sys::AHSSystem)
  G::Function = lrdf(sys)[1,1]
  ρ::Float64 = sys.ρ[1]

  F(s) = (s^2 * G(s) - 1) / s^3
  h(q) = -4π * real(F(im*q))

  S(q) = 1 + ρ*h(q)

  ret = Array{Function,2}(1,1)
  ret[1,1] = S

  return ret
end

function cavityfunction(sys::AHSSystem)
  @attach(sys, ρ, σ)

  # scalar constans
  ζ₂ = sum(ρ .* σ.^2)
  ζ₃ = sum(ρ .* σ.^3)
  η = (π/6) * ζ₃
  λ = 2π / (1-η)
  λ′ = π^2 * ζ₂ / (1-η)^2

  N = length(sys.σ)
  G = lrdf(sys)
  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    σᵢⱼ = (σ[i] + σ[j]) / 2

    g_raw = r -> inverselaplace(G[i,j], r, 128, 30) / r

    # Calibrate computational(?) error by an ugly way
    res = Optim.optimize(r -> -g_raw(r), 0.9σᵢⱼ, 1.1σᵢⱼ, rel_tol=1e-3)
    @assert Optim.converged(res)
    σ_raw = Optim.minimizer(res)

    Δσ = σ_raw - σᵢⱼ
    g(r) = r < σᵢⱼ ? 0. : g_raw(r + Δσ)

    ret[i,j] = spline(g, σᵢⱼ, R_MAX, 64, bc="extrapolate")
  end

  return ret
end

function prdf(ahs::AHSSystem)
  @attach(ahs, σ)

  N::Int = ncomp(ahs)
  y::Array{Function,2} = cavityfunction(ahs)

  ret = Array{Function}(N,N)

  for i in 1:N, j in 1:N
    σᵢⱼ = (σ[i] + σ[j]) / 2
    function g(r)::Float64
      if r < σᵢⱼ
        0
      else
        y[i,j](r)
      end
    end
    ret[i,j] = g
  end

  return ret
end

# Ref: K. Hoshino. and W. H. Young: J. Phys. F: Metal Phys. 10 (1980) 1365-1374.
function entropy(ahs::AHSSystem)::Float64
  @attach(ahs, σ)

  N::Int = ncomp(ahs)
  c::Vector{Float64} = composition(ahs)
  η::Float64= totalpackingfraction(ahs)
  ζ::Float64 = 1 / (1 - η)

  y₁::Float64 = 0.
  y₂::Float64 = 0.

  for i in 1:N, j in 1:N
    i ≥ j && continue

    σ_sum = sum(c .* σ)

    y₁ += c[i]*c[j] * (σ[i] + σ[j]) * (σ[i] - σ[j])^2 / σ_sum^3
    y₂ += c[i]*c[j] * σ[i]*σ[j] * (σ[i] - σ[j])^2 / σ_sum^6 * sum(c .* σ.^2)
  end

  #
  # Note that all S_xxx are divided by N*kB
  #
  S_η = - (ζ - 1)*(ζ + 3)

  S_σ = (3/2*(ζ^2 - 1) - 1/2*(ζ - 1)*(ζ - 3) - log(ζ))*y₁ + (3/2*(ζ - 1)^2 - 1/2*(ζ - 1)*(ζ - 3) - log(ζ))*y₂

  S = S_η + S_σ
end
