"""
nfe.jl

Types and functions for NFE (Nearly-Free Electrons) perturbing systems using
Bretonnet-Silbert psuedopotentials.

References:

* J. L. Bretonnet and A. Derouiche: Phys. Rev. B, 43 (1991), 8924-8929.
* N. Jakse and J. L. Bretonnet: J. Phys.: Condens. Matter, 7 (1995), 3803-3815.
"""

immutable NFESystem <: PertSystem
  ρ::Float64 # number density
  c::Vector{Float64} # composition
  T::Float64 # temperature
  m::Vector{Float64} # mass
  z::Vector{Float64} # number of electrons
  rc::Vector{Float64} # pseudopotential core radius
  a::Vector{Float64} # pseudopotential parameter
end

function NFESystem(ρ::Number, T::Number, m::Number, z::Number, rc::Number, a::Number)
  NFESystem(ρ, [1.0], T, [m], [z], [rc], [a])
end

function NFESystem(refsys::RefSystem; kwargs...)
end

function fermiwavenumber(sys::NFESystem)
  ρ = sys.ρ
  z̄ = sum(sys.c .* sys.z)

  return (3 * z̄ * ρ * π^2) ^ (1/3)
end

# form factor of pseudopotential (Bretonnet-Silbert)
function formfactor(sys::NFESystem)
  # @attach(sys, ρ, z, rc, a)
  ρ = sys.ρ

  N = length(sys.z)
  ω = Array{Function}(N)

  for i = 1:N
    z = sys.z[i]
    rc = sys.rc[i]
    a = sys.a[i]

    B₁ = (z/rc) * (1 - 2a/rc) * exp(rc/a)
    B₂ = (2z/rc) * (a/rc - 1) * exp(0.5rc/a)

    ω[i] = q -> begin
      J₁(q) = 2 - exp(-rc/a) * (
      (rc * (1 + a^2 * q^2) / a + 1 - a^2 * q^2) * sin(q*rc) / (a*q) +
      (2 + rc * (1 + a^2 * q^2) / a) * cos(q*rc) )

      J₂(q) = 2 - exp(-rc/2a) * (
      (rc * (1 + 4 * a^2 * q^2) / 2a + 1 - 4 * a^2 * q^2) * sin(q*rc) / (2a*q) +
      (2 + rc * (1 + 4 * a^2 * q^2) / 2a) * cos(q*rc) )

      4π * ρ * a^3 * (
      B₁ * J₁(q) / (1 + a^2 * q^2)^2 + 8B₂ * J₂(q) / (1 + 4 * a^2 * q^2)^2) -
      (4π * z * ρ / q^2) * cos(q*rc)
    end
  end

  return ω
end

# Hartree dielectric function
function dielectric(sys::NFESystem)
  ρ = sys.ρ
  z̄ = sum(sys.c .* sys.z)

  kF = (3 * z̄ * ρ * π^2) ^ (1/3) # Fermi wavenumber

  return ϵ(q) = begin
    x = q / 2kF
    1 + 1/(2π*kF) * (1/x^2) * (1 + (1-x^2)/2x * log(abs((1+x)/(1-x))))
  end
end

# local-field exchange-correlation function (Vashishta-Singwi)
function localfiled(sys::NFESystem)
  ρ = sys.ρ
  z̄ = sum(sys.c .* sys.z)

  kF = (3 * z̄ * ρ * π^2) ^ (1/3) # Fermi wavenumber
  rs = ((3/4π) / ρ / z̄)^(1/3) # electron distance

  A = 0.5362 + 0.1874rs - 0.0157rs^2 + 0.0008rs^3
  B = 0.42 - 0.0582rs + 0.0079rs^2 - 0.0005rs^3

  G(q) = A*(1 - exp(-B*(q/kF)^2))

  return G
end

# normalized wavenumber-energy characteristic
function wnechar(sys::NFESystem)
  @attach(sys, ρ, z)

  ω = formfactor(sys)
  ϵ = dielectric(sys)
  G = localfiled(sys)

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

function pairpotential(sys::NFESystem)
  @attach(sys, ρ, c, z, rc, a)

  z̄ = sum(c.*z)

  F = wnechar(sys)

  N = length(c)
  u = Array{Function}(N,N)

  for n in 1:length(u)
    i, j = ind2sub((N,N), n)
    u[i,j] = r -> begin
      𝐹(r) = ∫(q -> F[i,j](q) * sin(q*r)/q, eps(Float64), 20)
      z[i]*z[j] / r * (1 - 2/π * 𝐹(r))
    end
  end

  return u
end
