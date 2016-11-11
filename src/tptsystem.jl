type TPTSystem{Tr <: ReferenceSystem, Tp <: Perturbation}
  ncomp::Int # number of components
  density::Float64 # total number density
  composition::Vector{Float64}
  temp::Float64 # temperature
  mass::Vector{Float64} # mass of components
  ref::Tr # reference system
  pert::Tp # perturbing system
end

function TPTSystem(ref::IndependentReferenceSystem, pert::Perturbation = NullPerturbation(); T = 0, m = 0)
  N::Int = ncomp(ref)
  ρ::Float64 = totalnumberdensity(ref)
  c::Vector{Float64} = composition(ref)

  temp::Float64 = T
  mass::Vector{Float64} = []

  if typeof(m) <: Vector
    mass = m
  else m > 0
    mass = [m]
  end

  TPTSystem(N, ρ, c, temp, mass, ref, pert)
end

function prdf(sys::TPTSystem)
  prdf(sys.ref)
end

function psf(sys::TPTSystem)
  psf(sys.ref)
end

function pairpotential(sys::TPTSystem)
  pairpotential(sys.pert)
end

function freeenergy(sys::TPTSystem)
  if typeof(sys.ref) <: IndependentReferenceSystem
    F₀, rep0 = freeenergy(sys.ref)
  else
    F₀, rep0 = freeenergy(sys.ref, sys.pert)
  end
  F₁, rep1 = freeenergy(sys.ref, sys.pert)
  return F₀ + F₁, rep0, rep1
end

function potentialenergy(sys::TPTSystem)
  return potentialenergy(sys.ref, sys.pert)
end

function ncomp(sys::TPTSystem)::Int
  sys.ncomp
end

function totalnumberdensity(sys::TPTSystem)::Float64
  sys.density
end

function composition(sys::TPTSystem)::Vector{Float64}
  sys.composition
end

function temperature(sys::TPTSystem)::Float64
  sys.temp
end

function mass(sys::TPTSystem)::Vector{Float64}
  sys.mass
end

function kinetic(sys::TPTSystem)::Float64
  T = temperature(sys)
  K = 3/2 * kB*T
end

function entropy(sys::TPTSystem)::Float64
  ρ::Float64 = totalnumberdensity(sys)
  c::Vector{Float64} = composition(sys)
  T::Float64 = temperature(sys)
  m::Vector{Float64} = mass(sys)

  # Ideal gas entropy
  if m == Float64[] || T == 0.
    S_gas = 0
  else
    S_gas = 5/2 + log((prod(m.^c) * kB*T / 2π)^(3/2) / ρ)
  end

  S_conf = - sum(c .* log(c)) # configurational entropy

  S_ref = entropy(sys.ref) # reference system entropy
  S_pert = entropy(sys.pert, sys.ref) # perturbation entropy

  S = S_gas + S_conf + S_ref + S_pert
end

function helmholtz(sys::TPTSystem)::Float64
  T = temperature(sys)
  K = kinetic(sys)
  S = entropy(sys)

  F = K - kB*T*S
end
