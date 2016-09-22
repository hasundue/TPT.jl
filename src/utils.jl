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
