type TPTSystem{Tr <: ReferenceSystem, Tp <: Perturbation}
  ncomp::Int # number of components
  density::Float64 # total number density
  composition::Vector{Float64}
  temp::Float64 # temperature
  mass::Vector{Float64} # mass of components
  ref::Tr # reference system
  pert::Tp # perturbing system
end

function TPTSystem(ref::ReferenceSystem, pert::Perturbation = NullPerturbation(ncomp(ref)); T = InvalTemp, m = InvalMass)
  N::Int = ncomp(ref)
  ρ::Float64 = totalnumberdensity(ref)
  c::Vector{Float64} = composition(ref)

  temp::Float64 = T

  if temp == InvalTemp
    temp = temperature(ref) # may be also InvalTemp
  end

  mass::Vector{Float64} = []

  if typeof(m) <: Vector
    mass = m
  else m > 0
    mass = [m]
  end

  TPTSystem(N, ρ, c, temp, mass, ref, pert)
end

paircorrelation(sys::TPTSystem) = paircorrelation(sys.ref)
structurefactor(sys::TPTSystem) = structurefactor(sys.ref)
pairpotential(sys::TPTSystem) = pairpotential(sys.pert)
ncomp(sys::TPTSystem) = sys.ncomp
totalnumberdensity(sys::TPTSystem) = sys.density
composition(sys::TPTSystem) = sys.composition
temperature(sys::TPTSystem) = sys.temp
mass(sys::TPTSystem) = sys.mass

function kinetic(sys::TPTSystem)::Float64
  T = temperature(sys)
  K = 3/2 * kB*T
end

function entropy(sys::TPTSystem)::Float64
  T::Float64 = temperature(sys)

  S_gas = entropy_gas(sys)
  S_conf = entropy_conf(sys)
  S_ref = entropy(sys.ref) # reference system entropy
  S_pert = entropy(sys.pert, sys.ref, T) # perturbation entropy

  S = S_gas + S_conf + S_ref + S_pert
end

function entropy_gas(sys::TPTSystem)::Float64
  N::Int = ncomp(sys)
  ρ::Float64 = totalnumberdensity(sys)
  c::Vector{Float64} = composition(sys)
  T::Float64 = temperature(sys)
  m::Vector{Float64} = mass(sys)

  # Ideal gas entropy
  if m == InvalMass || T == InvalTemp
    return S = 0
  else
    P::Float64 = 1
    for i in 1:N
      if c[i] > 0
        P *= m[i]^c[i]
      end
    end
    return S = 5/2 + log((P * kB*T / 2π)^(3/2) / ρ)
  end
end

function entropy_conf(sys::TPTSystem)::Float64
  N::Int = ncomp(sys)
  c::Vector{Float64} = composition(sys)
  T::Float64 = temperature(sys)

  if T == InvalTemp
    return S = 0
  else
    S::Float64 = 0
    for i in 1:N
      if c[i] > 0
        S -= c[i]*log(c[i])
      end
    end
    return S
  end
end

function internal(sys::TPTSystem)::Float64
  U = internal(sys.pert, sys.ref)
end

function helmholtz(sys::TPTSystem)::Float64
  T = temperature(sys)
  K = kinetic(sys)
  U = internal(sys)
  S = entropy(sys)

  F = K + U - kB*T*S
end
