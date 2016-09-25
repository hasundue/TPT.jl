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
  a = 30.0
  h = π/t
  u = exp(a/2) / t
  x = a / 2t

  c(m::Int) = binomial(meuler, m-1)
  fs(s::Complex{Float64}) = real(F(s))

  suma = fs(x + 0im)/2
  for n in 1:nterm
    suma = suma + (-1)^n * fs(x + im*n*h)
  end

  su = Array{Float64}(meuler+2)
  su[1] = suma
  for k in 1:meuler+1
    n = nterm + k
    su[k+1] = su[k] + (-1)^n * fs(x + im*n*h)
  end

  argsu = 0
  argsu1 = 0
  for j in 1:meuler+1
    argsu = argsu + c(j)*su[j]
    argsu1 = argsu1 + c(j)*su[j+1]
  end

  fun = u * argsu / 2^meuler
  fun1 = u * argsu1 / 2^meuler

  return fun1
end

function tablize(f::Function, a::Float64, b::Float64, N::Int) :: Function
  @assert a < b "invalid arguments"

  Δx = (b - a) / N
  table = Vector{Float64}(N)
  for i in 1:N
    table[i] = f(a + i*Δx)
  end

  function ftab(x)
    if x < a || x > b
      error("the variable is out of the domain")
    end
    i = convert(Int, div(x - a, Δx))
    dx = rem(x - a, Δx)
    return table[i] + dx * (table[i+1] - table[i]) / Δx
  end

  return ftab
end
