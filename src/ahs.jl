immutable AHS <: RefSystem
  σ::Vector{Float64} # hard-sphere diameters
  ρ::Vector{Float64} # number densities
end # type

@doc doc"""
AHS(; kwargs...)

Initializes an additive hard-sphere reference system. You may use one of possible combination of keyword arguments.

One-component case:

* two out of σ, ρ, and η (all scalar)

Multi-component case:

* σ (vector) and ρ (vector)
* ρ (scalar), σ (vector), and c (vector)
* η (scalar), σ (vector), and c (vector)

You may use sigma, rho, or eta instead of σ, ρ, or η, respectively.
""" ->
function AHS(; kwargs...)
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

    return AHS([σ], [ρ])

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
    return AHS(σ, ρ)
  end
end

# ref: S. B. Yuste et al: J. Chem. Phys., Vol. 108, No.9, 1 (1998), 3683-3693.
function lrdf(sys::AHS)
  if length(sys.σ) == 1
    lrdf1(sys)
  elseif length(sys.σ) == 2
    lrdf2(sys)
  else
    error("unsupported number of components")
  end
end

function lrdf1(sys::AHS)
  σ = sys.σ[1]
  ρ = sys.ρ[1]
  η = π/6 * ρ * σ^3

  if η > π/(3*√2)
    warn("unphysical packing fraction (η > 0.74)")
  end

  L⁰ = 2π * (1 + 2η) / (1 - η)^2
  L¹ = 2π*σ * (1 + 1/2*η) / (1 - η)^2

  ϕ₁(x) = x^-2 * (1 - x - exp(-x))
  ϕ₂(x) = x^-3 * (1 - x + x^2/2 - exp(-x))

  G(s) = exp(-σ*s) / (2π*s^2) * (L⁰ + L¹*s) /
         (1 - ρ*(ϕ₂(σ*s)*σ^3*L⁰ + ϕ₁(σ*s)*σ^2*L¹))

  return G
end

function lrdf2(sys::AHS)
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
psf(system::AHS)

Returns a structure factor of additive hard-sphere system in wave-number space. Currently the cases where N = 1 or N = 2 are only supported.

ref: S. B. Yuste et al: J. Chem. Phys., Vol. 108, No.9, 1 (1998), 3683-3693.
""" ->
function psf(sys::AHS)
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

function psf1(sys::AHS)
  G = lrdf(sys)
  ρ = sys.ρ[1]

  F(s) = (s^2 * G(s) - 1) / s^3
  h(q) = -4π * real(F(im*q))

  S(q) = 1 + ρ*h(q)

  return S
end
