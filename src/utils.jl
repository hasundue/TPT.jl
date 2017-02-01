function expandkwargs(kwargs::Vector)
  m = length(kwargs)
  keys = Vector{Symbol}(m)
  vals = Vector{Any}(m)
  for i in 1:m
    keys[i] = kwargs[i][1]
    vals[i] = kwargs[i][2]
  end
  return keys, vals
end

macro attach(obj, fields...)
  block = Expr(:block)
  for f in fields
    e = :($f = $obj.$f)
    push!(block.args, e)
  end
  return esc(:($block))
end

const greekdict = Dict(:σ => :sigma, :ρ => :rho, :η => :eta)

function greekin(g::Symbol, c::AbstractArray)
  if g in c
    return true
  else
    s = get(greekdict, g, :notfound)
    if s == :notfound
      error("an unknown greek letter is given")
    else
      return s in c
    end
  end
end

function greekfind(c::AbstractArray, g::Symbol)
  if g in c
    return findfirst(c, g)
  else
    s = get(greekdict, g, :notfound)
    if s == :notfound
      error("an unknown greek letter is given")
    else
      return findfirst(c, s)
    end
  end
end

function ∫(f::Function, a::Number, b::Number; e::Float64=1e-6)
  val, err = quadgk(f, a, b, reltol=e)
  return val
end

# Ref: http://www.eweb.unex.es/eweb/fisteor/santos/files/rdf.nb
function inverselaplace(F::Function, t::Float64, nterm::Int, meuler::Int)
  a::Float64 = 30.0
  h::Float64 = π/t
  u::Float64 = exp(a/2) / t
  x::Float64 = a / 2t

  c(m::Int)::Int = binomial(meuler, m-1)
  fs(s::Complex{Float64})::Float64 = real(F(s))

  suma::Float64 = fs(x + 0im)/2
  @simd for n in 1:nterm
    @inbounds suma += (-1)^n * fs(x + im*n*h)
  end

  su = Vector{Float64}(meuler+2)
  su[1] = suma
  for k in 1:meuler+1
    n = nterm + k
    @inbounds su[k+1] = su[k] + (-1)^n * fs(x + im*n*h)
  end

  # argsu = 0
  argsu1::Float64 = 0
  @simd for j in 1:meuler+1
    # argsu = argsu + c(j)*su[j]
    @inbounds argsu1 += c(j)*su[j+1]
  end

  # fun = u * argsu / 2^meuler
  fun1::Float64 = u * argsu1 / 2^meuler

  return fun1
end

function spline(f::Function, a::Float64, b::Float64, N::Int; bc="error")::Function
  spl = spline_spl(f, a, b, N, bc)

  fspl(x::Real)::Float64 = evaluate(spl, x)
end

function spline_derivative(f::Function, a::Float64, b::Float64, N::Int;
                           bc="error")::Function
  spl = spline_spl(f, a, b, N, bc)

  f′spl(x::Real)::Float64 = derivative(spl, x)
end

function spline_integral(f::Function, a::Float64, b::Float64, N::Int)::Float64
  spl = spline_spl(f, a, b, N, "zero")
  integrate(spl, a, b)
end

function spline_spl(f::Function, a::Float64, b::Float64, N::Int,
                    bc::AbstractString)::Spline1D
  @assert a < b "invalid arguments for spline"

  Δx::Float64 = (b - a) / N
  𝐱::Vector{Float64} = collect(a : Δx : b)
  𝐲::Vector{Float64} = [f(x) for x in 𝐱]

  Spline1D(𝐱, 𝐲, k=3, bc=bc)
end

function cardano(a₃::Real, a₂::Real, a₁::Real, a₀::Real)
  A₂ = a₂ / a₃
  A₁ = a₁ / a₃
  A₀ = a₀ / a₃

  p = (A₁ - 1/3*A₂^2) / 3
  q = (A₀ - 1/3*A₁*A₂ + 2/27*A₂^3) / 2

  D = - (q^2 + p^3)

  if D ≈ 0 && q ≈ 0
    y = [0]
  elseif D ≈ 0
    y = [-2*∛q, ∛q]
  elseif D > 0
    Θ = atan2(√D, -q)
    R = (q^2 + D)^(1/6) * cos(Θ/3)
    Q = (q^2 + D)^(1/6) * sin(Θ/3)
    y = [2R, -R - √3*Q, -R + √3*Q]
  elseif D < 0
    S = ∛(-q + √(-D))
    T = ∛(-q - √(-D))
    y = [S+T, -1/2*(S+T) + √3/2*im*(S-T), -1/2*(S+T) - √3/2*im*(S-T)]
  end

  x = y - A₂/3
end
