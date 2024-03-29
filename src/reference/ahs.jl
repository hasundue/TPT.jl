"""
ahs.jl

Types and functions for additive hard-sphere reference systems.

References:

* S. B. Yuste et al: J. Chem. Phys., Vol. 108, No.9, 1 (1998), 3683-3693.
* K. Hoshino. and W. H. Young: J. Phys. F: Metal Phys., Vol. 10 (1980), 1365-1374.

"""

immutable AHS <: IndependentReferenceSystem
  approx::Symbol # approximation
  σ::Vector{Float64} # hard-sphere diameters
  ρ::Vector{Float64} # number densities
  α::Float64 # α parameter for RFA approximation
end # type

function AHS(approx::Symbol, σ::Vector{Float64}, ρ::Vector{Float64})
  optimizealpha(AHS(approx, σ, ρ, 0.))
end

@doc doc"""
AHS(; kwargs...)

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
function AHS(; kwargs...)
  keys, vals = expandkwargs(kwargs)

  # The type of approximation
  if in(:approx, keys)
    ind = findfirst(keys, :approx)
    approx = vals[ind]

    @assert in(approx, [:PY, :RFA]) "unsuportted type of approximation for AHS"

    deleteat!(keys, ind)
    deleteat!(vals, ind)
  else
    approx = :RFA # default
  end

  @assert 2 ≤ length(keys) ≤ 3 "wrong number of arguments for AHS"

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
      error("invalid arguments for AHS")
    end # if

    return AHS(approx, [σ], [ρ])
  else # Multi-component case

    if length(kwargs) == 2
      @assert all(x-> typeof(x) <: Vector, vals) "invalid arguments for AHS"
      σ = vals[greekfind(keys, :σ)]
      ρ = vals[greekfind(keys, :ρ)]

    else
      @assert greekin(:σ, keys) && greekin(:c, keys) "invalid arguments for AHS"
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

    @assert length(σ) == length(ρ) "invalid arguments for AHS"

    return AHS(approx, σ, ρ)
  end
end

ncomp(ahs::AHS)::Int = length(ahs.σ)
composition(ahs::AHS)::Vector{Float64} = ahs.ρ / sum(ahs.ρ)
numberdensity(ahs::AHS)::Vector{Float64} = ahs.ρ
totalnumberdensity(ahs::AHS)::Float64 = sum(ahs.ρ)
packingfraction(ahs::AHS)::Vector{Float64} = (π/6) * ahs.ρ .* ahs.σ.^3
totalpackingfraction(ahs::AHS)::Float64 = sum(packingfraction(ahs))
temperature(ahs::AHS)::Float64 = InvalTemp

function hsdiameter(sys::AHS)::Array{Float64,2}
  N = ncomp(sys)

  ret = Array{Float64,2}(N,N)

  for i in 1:N, j in 1:N
    if i == j
      ret[i,j] = sys.σ[i]
    else
      ret[i,j] = (sys.σ[i] + sys.σ[j]) / 2
    end
  end

  return ret
end

emptyradius(ahs::AHS)::Array{Float64,2} = hsdiameter(ahs)

function pairpotential(ahs::AHS)::Array{Function,2}
  N = ncomp(ahs)
  σ = hsdiameter(ahs)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    function u(r)::Float64
      if r < σ[i,j]
        Inf
      else
        0
      end
    end

    ret[i,j] = ret[j,i] = u
  end

  return ret
end

function zeta(ahs::AHS, n::Int)::Float64
  N = ncomp(ahs)
  sum(ahs.ρ[i] * ahs.σ[i]^n for i = 1:N)
end

@doc doc"""
pressure(ahs::AHS)::Float64

Returns pressure of an additive hard-sphere system multiplied by inverse
temperature based on BMCSL equation of state.

""" ->
function pressure(ahs::AHS)::Float64
  ρ = totalnumberdensity(ahs)
  η = totalpackingfraction(ahs)
  ζ₁ = zeta(ahs, 1)
  ζ₂ = zeta(ahs, 2)

  ρ/(1-η) + π/2*ζ₁*ζ₂/(1-η)^2 + π^2/36*ζ₂^3*(3-η)/(1-η)^3
end

@doc doc"""
chemicalpotential(ahs::AHS)::Vector{Float64}

Returns a vector of excess (residual) chemical potential of each specie
multiplied by inverse temperature.

""" ->
function chemicalpotential(ahs::AHS)::Vector{Float64}
  N = ncomp(ahs)
  c = composition(ahs)
  σ = ahs.σ
  ρ = totalnumberdensity(ahs)

  R = σ/2
  S = π * σ.^2
  V = π/6 * σ.^3

  v = ρ * sum(c .* V)
  s = ρ * sum(c .* S)
  q = ρ * sum(c .* R.^2)
  r = ρ * sum(c .* R)

  βΔμ = Vector{Float64}(N)
  for i in 1:N
    βΔμ[i] = -log(1-v) + (R[i]*s + S[i]*r + V[i]*ρ)/(1-v) +
             (R[i]^2*s^2 + 2S[i]*q*s + 6V[i]*r*s)/6/(1-v)^2 +
             V[i]*q*s^2*(2 - v/3)/9/(1-v)^3
  end

  return βΔμ
end

function chemicalpotentialPYC(ahs::AHS)::Vector{Float64}
  N = ncomp(ahs)
  σ = hsdiameter(ahs)
  ρ = numberdensity(ahs)
  η = π/6*ρ
  ξ = sum( η[i]*σ[i,i]^3 for i in 1:N )
  X = sum( η[i]*σ[i,i]^2 for i in 1:N )
  Y = sum( η[i]*σ[i,i]^1 for i in 1:N )

  βp = ( sum(ρ)*(1 + ξ + ξ^2) -
         18/π * sum( η[i]*η[j]*(σ[i,i] - σ[j,j])^2*(2σ[i,j] + σ[i,i]*σ[j,j]*X)
                     for i in 1:N, j in 1:N if i < j )
       ) / (1-ξ)^3

  βΔμ = Vector{Float64}(N)
  for i in 1:N
    βΔμ[i] = log(1-ξ) + π/6*βp*σ[i,i]^3 + σ[i,i]^2/(1-ξ)^3 *
             ( 3Y*(1-2ξ) + 9/2*X^2*(1-ξ) + 3ξ^2*Y ) +
             σ[i,i]/(1-ξ)^3*(3X - 6ξ*X + 3ξ^2*X)
  end

  return βΔμ
end

# Ref: A. Meyer et al. Chemical Physics 49 (1980) 147-152
function contactvaluePY(ahs::AHS)::Matrix{Float64}
  N = ncomp(ahs)
  if N > 1
    error("calculation of contact values of multi-component system in PY approxmation is not implemented")
  end
  η = totalpackingfraction(ahs)
  ξ₀ = (1 + η/2) / (1 - η)^2
  g₀ = Matrix{Float64}(N,N)
  g₀[1,1] = ξ₀
  g₀
end

function contactvalueRFA(ahs::AHS)::Matrix{Float64}
  N::Int = ncomp(ahs)
  σ::Array{Float64,2} = hsdiameter(ahs)
  ρ::Vector{Float64} = numberdensity(ahs)

  # scalar constans
  ζ₂::Float64 = sum(ρ[i] * σ[i,i]^2 for i = 1:N)
  ζ₃::Float64 = sum(ρ[i] * σ[i,i]^3 for i = 1:N)
  η::Float64 = (π/6) * ζ₃
  λ::Float64 = 2π / (1-η)
  λ′::Float64 = π^2 * ζ₂ / (1-η)^2

  [ 1/2π * (λ + λ′/2 * σ[i,i]*σ[j,j]/σ[i,j] +
            λ′^2 / 18λ * σ[i,i]^2 * σ[j,j]^2 / σ[i,j]^2)
    for i = 1:N, j = 1:N ]
end

function contactvalue(ahs::AHS)::Matrix{Float64}
  if ahs.approx == :PY
    contactvaluePY(ahs)
  elseif ahs.approx == :RFA
    contactvalueRFA(ahs)
  end
end

function contactgradientPY(ahs::AHS)::Matrix{Float64}
  N = ncomp(ahs)
  if N > 1
    error("calculation of contact gradients of multi-component system in PY approxmation is not implemented")
  end
  η = totalpackingfraction(ahs)
  ξ₀ = (1 + η/2) / (1 - η)^2
  ξ₁ = 9/2 * η * (1 + η) / (1 - η)^3
  Y = Matrix{Float64}(N,N)
  Y[1,1] = ξ₀ / ξ₁
  Y
end

function contactgradientRFA(ahs::AHS)::Matrix{Float64}
  N::Int = ncomp(ahs)
  σ::Array{Float64,2} = hsdiameter(ahs)
  ρ::Vector{Float64} = numberdensity(ahs)
  α::Float64 = ahs.α
  g_σ::Array{Float64,2} = contactvalue(ahs)

  ζ₂::Float64 = sum(ρ[i] * σ[i,i]^2 for i = 1:N)
  ζ₃::Float64 = sum(ρ[i] * σ[i,i]^3 for i = 1:N)
  η::Float64 = (π/6) * ζ₃
  λ::Float64 = 2π / (1-η)
  λ′::Float64 = π^2 * ζ₂ / (1-η)^2

  L²::Array{Float64,2} =
    [ 2π * α * σ[i,j] * g_σ[i,j] for i = 1:N, j = 1:N ]

  L¹::Array{Float64,2} =
    [ ( λ*σ[i,j] + λ′/2*σ[i,i]*σ[j,j] + (λ + λ′*σ[i,i])*α -
        λ/2*σ[i,i]*sum(ρ[k]*σ[k,k]*L²[k,j] for k = 1:N) )
      for i = 1:N, j = 1:N ]

  g′::Array{Float64,2} =
    [ 1 / (2π*α*σ[i,j]) * (L¹[i,j] - L²[i,j]*(1/α + 1/σ[i,j]))
      for i in 1:N, j in 1:N ]
end

function contactgradient(ahs::AHS)::Matrix{Float64}
  if ahs.approx == :PY
    contactgradientPY(ahs)
  elseif ahs.approx == :RFA
    contactgradientRFA(ahs)
  end
end

#
# Pair correlation function in Laplace space
# ref: S. B. Yuste et al: J. Chem. Phys., Vol. 108, No.9, 1 (1998), 3683-3693.
#
function optimizealpha(ahs::AHS)
  ahs.approx == :PY && return ahs

  N::Int = ncomp(ahs)
  x::Vector{Float64} = composition(ahs)
  σ::Array{Float64,2} = hsdiameter(ahs)
  ρ::Vector{Float64} = numberdensity(ahs)

  # scalar constans
  ζ₁::Float64 = sum(ρ[i] * σ[i,i] for i = 1:N)
  ζ₂::Float64 = sum(ρ[i] * σ[i,i]^2 for i = 1:N)
  ζ₃::Float64 = sum(ρ[i] * σ[i,i]^3 for i = 1:N)
  η::Float64 = (π/6) * ζ₃
  λ::Float64 = 2π / (1-η)
  λ′::Float64 = π^2 * ζ₂ / (1-η)^2

  # auxiliary qunatities and functions
  I::Array{Float64,2} = δ::Array{Float64,2} = eye(N)
  ϕ₀(x::Complex{Float64}) = x^-1 * (1 - exp(-x))
  ϕ₁(x::Complex{Float64}) = x^-2 * (1 - x - exp(-x))
  ϕ₂(x::Complex{Float64}) = x^-3 * (1 - x + x^2/2 - exp(-x))

  # isothermal susceptibity
  χ::Float64 = sum(ρ) / (sum(ρ)/(1-η)^2 + π*ζ₁*ζ₂/(1-η)^3 +
                         π^2/36 * ζ₂^3 * (9-4η+η^2)/(1-η)^4)

  # the contact value of pair correlation function
  g_σ::Array{Float64,2} = contactvalue(ahs)

  L⁰ = Array{Float64,2}(N,N)
  L¹ = Array{Float64,2}(N,N)
  L² = Array{Float64,2}(N,N)
  A⁰ = Array{Float64,2}(N,N)
  A¹ = Array{Float64,2}(N,N)
  A² = Array{Float64,2}(N,N)
  A³ = Array{Float64,2}(N,N)
  B⁰ = Array{Float64,2}(N,N)
  B¹ = Array{Float64,2}(N,N)
  H⁰ = Array{Float64,2}(N,N)
  H¹ = Array{Float64,2}(N,N)

  function fopt(v::Vector{Float64}, g::Vector)::Float64
    α = v[1]

    # express L² in terms of α
    L² = [ 2π * α * σ[i,j] * g_σ[i,j] for i = 1:N, j = 1:N ]

    #
    # express H¹ in terms of L² and α
    #
    for i in 1:N, j in 1:N
      L⁰[i,j] = λ + λ′*σ[j,j] + 2λ′*α - λ*sum(ρ[k]*σ[k,k]*L²[k,j] for k = 1:N)
      L¹[i,j] = λ*σ[i,j] + λ′/2*σ[i,i]*σ[j,j] + (λ + λ′*σ[i,i])*α -
                λ/2*σ[i,i]*sum(ρ[k]*σ[k,k]*L²[k,j] for k = 1:N)

      Aᵢⱼ(n)::Float64 = (-1)^n * ρ[i] *
        ( σ[i,i]^(n+3) / factorial(n+3) * L⁰[i,j] -
          σ[i,i]^(n+2) / factorial(n+2) * L¹[i,j] +
          σ[i,i]^(n+1) / factorial(n+1) * L²[i,j] )

      A⁰[i,j] = Aᵢⱼ(0)
      A¹[i,j] = Aᵢⱼ(1)
      A²[i,j] = Aᵢⱼ(2)
      A³[i,j] = Aᵢⱼ(3)
    end

    for i in 1:N, j in 1:N
      B⁰[i,j] = 1/2π * L²[i,j] + sum(A²[k,j] for k = 1:N) -
        sum(σ[i,k]*(α*δ[k,j] - A¹[k,j]) for k = 1:N) -
        sum(1/2 * σ[i,k]^2 * (δ[k,j] - A⁰[k,j]) for k = 1:N)
    end

    H⁰ = B⁰ * inv(I - A⁰)

    for i in 1:N, j in 1:N
      B¹[i,j] = sum(A³[k,j] for k = 1:N) + sum(σ[i,k]*A²[k,j] for k = 1:N) -
        sum((1/2*σ[i,k]^2 + H⁰[i,k])*(α*δ[k,j] - A¹[k,j]) for k = 1:N) -
        sum((1/6*σ[i,k]^3 + σ[i,k]*H⁰[i,k])*(δ[k,j] - A⁰[k,j]) for k = 1:N)
    end

    H¹ = B¹ * inv(I - A⁰)

    #
    # calculate inversed isothermal susceptibity for given α
    #
    h0::Array{Float64,2} =
      [ -4π * √(ρ[i]*ρ[j]) * H¹[i,j] for i = 1:N, j = 1:N ]

    χ⁻¹ = sum(√(x[i]*x[j]) * inv(I + h0)[i,j] for i = 1:N, j = 1:N)

    abs(χ - 1/χ⁻¹)
  end # fopt

  # heuristic initial and boundary condition
  σ_m::Float64 = sum(x .* σ)
  α_init::Float64 = 0.01 * σ_m / η
  α_min::Float64 = 0.0 * α_init
  α_max::Float64 = 3.0 * α_init

  # opt = Optim.optimize(fopt, α_min, α_max, abs_tol=1e-6)
  # show(opt)
  # α::Float64 = Optim.minimizer(opt)

  opt = NLopt.Opt(:LN_BOBYQA, 1)
  NLopt.min_objective!(opt, fopt)
  NLopt.xtol_rel!(opt, 1e-4)
  NLopt.lower_bounds!(opt, [α_min])
  NLopt.upper_bounds!(opt, [α_max])

  (fmin, xmin, res) = NLopt.optimize(opt, [α_init])
  α::Float64 = xmin[1]

  if isapprox(α, α_min, rtol=1e-2) || isapprox(α, α_max, rtol=1e-2)
    warn("α = $α is on a bound of optimization (α_init = $α_init)")
  end

  AHS(ahs.approx, ahs.σ, ahs.ρ, α)
end # optimizealpha()

#
# Pair correlation function in Laplace space
# ref: S. B. Yuste et al: J. Chem. Phys., Vol. 108, No.9, 1 (1998), 3683-3693.
#
function paircorrelationlaplace(ahs::AHS)::Array{Function,2}
  N::Int = ncomp(ahs)
  x::Vector{Float64} = composition(ahs)
  ρ::Vector{Float64} = numberdensity(ahs)
  σ::Array{Float64,2} = hsdiameter(ahs)
  g_σ::Array{Float64,2} = contactvalue(ahs)

  # scalar constans
  ζ₂::Float64 = sum(ρ[i] * σ[i,i]^2 for i = 1:N)
  ζ₃::Float64 = sum(ρ[i] * σ[i,i]^3 for i = 1:N)
  η::Float64 = (π/6) * ζ₃
  λ::Float64 = 2π / (1-η)
  λ′::Float64 = π^2 * ζ₂ / (1-η)^2

  # auxiliary qunatities and functions
  I::Array{Float64,2} = δ::Array{Float64,2} = eye(N)
  ϕ₀(x::Complex{Float64}) = x^-1 * (1 - exp(-x))
  ϕ₁(x::Complex{Float64}) = x^-2 * (1 - x - exp(-x))
  ϕ₂(x::Complex{Float64}) = x^-3 * (1 - x + x^2/2 - exp(-x))

  α::Float64 = ahs.α

  L²::Array{Float64,2} =
    [ 2π * α * σ[i,j] * g_σ[i,j] for i = 1:N, j = 1:N ]

  L⁰::Array{Float64,2} =
    [ ( λ + λ′*σ[j,j] + 2λ′*α - λ*sum(ρ[k]*σ[k,k]*L²[k,j] for k = 1:N) )
      for i = 1:N, j = 1:N ]

  L¹::Array{Float64,2} =
    [ ( λ*σ[i,j] + λ′/2*σ[i,i]*σ[j,j] + (λ + λ′*σ[i,i])*α -
        λ/2*σ[i,i]*sum(ρ[k]*σ[k,k]*L²[k,j] for k = 1:N) )
      for i = 1:N, j = 1:N ]

  𝐋::Array{Function,2} =
    [ s::Complex{Float64} -> L⁰[i,j] + L¹[i,j]*s + L²[i,j]*s^2
      for i = 1:N, j = 1:N ]

  𝐀::Array{Function,2} =
    [ s::Complex{Float64} -> ρ[i] * ( ϕ₂(σ[i,i]*s)*σ[i,i]^3*L⁰[i,j] +
                                      ϕ₁(σ[i,i]*s)*σ[i,i]^2*L¹[i,j] +
                                      ϕ₀(σ[i,i]*s)*σ[i,i]^1*L²[i,j] )
      for i = 1:N, j = 1:N ]

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    function G(s::Complex{Float64})::Complex{Float64}
      L::Array{Complex{Float64},2} = [ 𝐋[i,j](s) for i = 1:N, j = 1:N ]
      A::Array{Complex{Float64},2} = [ 𝐀[i,j](s) for i = 1:N, j = 1:N ]

      exp(-σ[i,j]*s) / (2π*s^2) * (L * inv((1 + α*s)*I - A))[i,j]
    end

    ret[i,j] = ret[j,i] = G
  end

  return ret
end


@doc doc"""
structurefactor(ahs::AHS)

Returns a structure factor of additive hard-sphere system in wave-number space.
""" ->
function structurefactor(ahs::AHS)::Array{Function,2}
  N::Int = ncomp(ahs)
  ρ::Vector{Float64} = numberdensity(ahs)

  # g(r)*r in Laplace space
  G::Array{Function,2} = paircorrelationlaplace(ahs)

  I = eye(N) # identity matrix

  # matrix of structure factors
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    F(s::Complex{Float64})::Complex{Float64} = (s^2 * G[i,j](s) - 1) / s^3
    hᵢⱼ(q::Float64)::Float64 = -4π * real(F(im*q))

    S(q::Real)::Float64 = I[i,j] + √(ρ[i]*ρ[j]) * hᵢⱼ(convert(Float64, q))

    ret[i,j] = ret[j,i] = S
  end

  return ret
end

function paircorrelation(ahs::AHS)::Array{Function,2}
  N::Int = ncomp(ahs)
  σ::Array{Float64,2} = hsdiameter(ahs)
  ρ::Vector{Float64} = numberdensity(ahs)

  # g(r)*r in Laplace space
  G::Array{Function,2} = paircorrelationlaplace(ahs)

  # g(r) at r = σ
  g_σ::Array{Float64,2} = contactvalue(ahs)

  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    function g(r)::Float64
      if r < σ[i,j]
        0
      elseif r == σ[i,j]
        g_σ[i,j]
      else
        nterm = trunc(Int, 2048/sqrt(σ[i,j]))
        inverselaplace(G[i,j], r, nterm, 30) / r
      end
    end

    ret[i,j] = ret[j,i] = g
  end

  return ret
end

function cavityfunction(ahs::AHS)
  N = ncomp(ahs)
  δ = eye(N,N)
  σ = hsdiameter(ahs)
  ρ = numberdensity(ahs)
  g₀ = contactvalue(ahs)

  M = chemicalpotential(ahs)
  S = [ -π*δ[i,j] * sum(ρ[k]*σ[k,j]^2*g₀[k,j] for k in 1:N) for i in 1:N, j in 1:N ]
  Y = contactvalue(ahs)
  Y′ = contactgradient(ahs)
  a = [ abs(σ[i,i] - σ[j,j]) / 2 for i in 1:N, j in 1:N ]

  ret = Matrix{Function}(N,N)
  for i in 1:N, j in 1:N
    i > j && continue
    s = σ[i,j] - a[i,j]
    α = Vector{Float64}(4)
    α[1] = M[i]
    α[2] = S[i,j]
    α[3] = s^-2 * ( 3log(Y[i,j]) - s*Y′[i,j]/Y[i,j] - 3α[1] - 2s*α[2] )
    α[4] = s^-3 * ( -2log(Y[i,j]) + s*Y′[i,j]/Y[i,j] + 2α[1] + s*α[2] )
    function y(r::Real)::Float64
      if r < a[i,j]
        exp(M[i])
      elseif r < σ[i,j]
        exp(sum( α[k]*(r - a[i,j])^(k-1) for k in 1:4 ))
      else
        paircorrelation(ahs)[i,j](r)
      end
    end
    ret[i,j] = ret[j,i] = y
  end

  return ret
end

# Ref: K. Hoshino. and W. H. Young: J. Phys. F: Metal Phys. 10 (1980) 1365-1374.
function entropy(ahs::AHS)::Float64
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

  S_η = - (ζ - 1)*(ζ + 3)

  S_σ = (3/2*(ζ^2 - 1) - 1/2*(ζ - 1)*(ζ - 3) - log(ζ))*y₁ + (3/2*(ζ - 1)^2 - 1/2*(ζ - 1)*(ζ - 3) - log(ζ))*y₂

  S = S_η + S_σ
end

# Perturbation-independent part of internal energy
function internal(ahs::AHS)::Float64
  U = 0
end
