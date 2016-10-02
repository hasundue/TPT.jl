"""
nfe.jl

Types and functions for NFE (Nearly-Free Electrons) perturbation.

References:

* J. L. Bretonnet and A. Derouiche: Phys. Rev. B, 43 (1991), 8924-8929.
* N. Jakse and J. L. Bretonnet: J. Phys.: Condens. Matter, 7 (1995), 3803-3815.
"""

const qmax = 20.0

immutable NFE{T <: PseudoPotential} <: NFEPerturbation
  ρ::Float64 # number density
  c::Vector{Float64} # composition
  T::Float64 # temperature
  m::Vector{Float64} # mass
  z::Vector{Float64} # number of electrons
  pp::T # pseudopotential
end

function NFE(ρ::Number, T::Number, m::Number, z::Number, pp::PseudoPotential)
  NFE(ρ, [1.0], T, [m], [z], pp)
end

function fermiwavenumber(nfe::NFE)
  ρ = nfe.ρ
  z̄ = sum(nfe.c .* nfe.z)

  return (3 * z̄ * ρ * π^2) ^ (1/3)
end

# Hartree dielectric function
function dielectric(nfe::NFE)
  ρ = nfe.ρ
  z̄ = sum(nfe.c .* nfe.z)

  kF = (3 * z̄ * ρ * π^2) ^ (1/3) # Fermi wavenumber

  return ϵ(q) = begin
    x = q / 2kF
    1 + 1/(2π*kF) * (1/x^2) * (1 + (1-x^2)/2x * log(abs((1+x)/(1-x))))
  end
end

# local-field exchange-correlation function (Vashishta-Singwi)
function localfiled(nfe::NFE)
  ρ = nfe.ρ
  z̄ = sum(nfe.c .* nfe.z)

  kF = (3 * z̄ * ρ * π^2) ^ (1/3) # Fermi wavenumber
  rs = ((3/4π) / ρ / z̄)^(1/3) # electron distance

  A = 0.5362 + 0.1874rs - 0.0157rs^2 + 0.0008rs^3
  B = 0.42 - 0.0582rs + 0.0079rs^2 - 0.0005rs^3

  G(q) = A*(1 - exp(-B*(q/kF)^2))

  return G
end

# normalized wavenumber-energy characteristic
function wnechar(nfe::NFE)
  @attach(nfe, ρ, z)

  ω = formfactor(nfe)
  ϵ = dielectric(nfe)
  G = localfiled(nfe)

  N = length(z)
  F = Array{Function}(N,N)

  for n in 1:length(F)
    i, j = ind2sub((N,N), n)
    F[i,j] = q -> begin
      q^4/(4π*ρ)^2 / (z[i]*z[j]) * (1 - 1/ϵ(q)) / (1 - G(q)) * ω[i](q) * ω[j](q)
    end
  end

  return F
end

function pairpotential(nfe::NFE)
  @attach(nfe, ρ, c, z)

  z̄ = sum(c.*z)

  F = wnechar(nfe)

  N = length(c)
  u = Array{Function}(N,N)

  for n in 1:length(u)
    i, j = ind2sub((N,N), n)
    u[i,j] = r -> begin
      𝐹(r) = ∫(q -> F[i,j](q) * sin(q*r)/q, eps(Float64), qmax)
      z[i]*z[j] / r * (1 - 2/π * 𝐹(r))
    end
  end

  return u
end
