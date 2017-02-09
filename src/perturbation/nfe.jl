"""
nfe.jl

Types and functions for NFE (Nearly-Free Electrons) perturbation.

References:

* J. L. Bretonnet and A. Derouiche: Phys. Rev. B, 43 (1991), 8924-8929.
* N. Jakse and J. L. Bretonnet: J. Phys.: Condens. Matter, 7 (1995), 3803-3815.
"""

immutable NFE{T <: PseudoPotential} <: NFEPerturbation
  ρ::Float64 # number density
  z::Float64 # mean number of valence electrons
  pseudo::T # pseudopotential
  localfield::Symbol # kind of local field function
end

function NFE(ρ::Real, z::Real, pseudo::PseudoPotential, lf::Symbol = :VS)
  NFE(ρ, z, pseudo, lf)
end

function NFE(ref::ReferenceSystem, pseudo::PseudoPotential, lf::Symbol = :VS)
  ρ::Float64 = totalnumberdensity(ref)

  c::Vector{Float64} = composition(ref)
  z::Float64 = sum(c .* pseudo.z)

  NFE(ρ, z, pseudo, lf)
end

function TPTSystem(ref::ReferenceSystem, pseudo::PseudoPotential; kwargs...)
  nfe = NFE(ref, pseudo)
  TPTSystem(ref, nfe; kwargs...)
end

ncomp(nfe::NFE)::Int = length(nfe.pseudo.rc)

function fermiwavenumber(nfe::NFE)::Float64
  kF = (3π^2 * nfe.z * nfe.ρ)^(1/3)
end

#
# Lindhard (Hartree) dielectric function
# Ref: W. A. Harrison: Elementary Electronic Structure Revised Edition, 462
#
function dielectric(nfe::NFE)::Function
  ρ::Float64 = nfe.ρ

  kF::Float64 = fermiwavenumber(nfe)

  ϵ(q)::Float64 = begin
    x = q / 2kF
    1 + 1/(2π*kF) * (1/x^2) * (1 + (1-x^2)/2x * log(abs((1+x)/(1-x))))
  end
end

function screenedformfactor(nfe::NFE)::Vector{Function}
  N::Int = length(nfe.pseudo.z)
  ω₀::Vector{Function} = formfactor(nfe)
  ϵ::Function = dielectric(nfe)

  ω = Vector{Function}(N)

  for i in 1:N
    ω[i] = q -> ω₀[i](q) / ϵ(q)
  end

  ω
end

function electrondistance(nfe::NFE)::Float64
  @attach(nfe, ρ, z)
  rs = ((3/4π) / ρ / z)^(1/3)
end

function localfield(nfe::NFE)::Function
  if nfe.localfield == :IU
    localfield_IU(nfe)
  elseif nfe.localfield == :VS
    localfield_VS(nfe)
  else
    error("NFE: unsupported type of local field function")
  end
end

# local-field exchange-correlation function (Vashishta-Singwi)
function localfield_VS(nfe::NFE)::Function
  ρ::Float64 = nfe.ρ
  z̄::Float64 = nfe.z

  kF::Float64 = fermiwavenumber(nfe)
  rs::Float64 = electrondistance(nfe)

  A::Float64 = 0.5362 + 0.1874rs - 0.0157rs^2 + 0.0008rs^3
  B::Float64 = 0.42 - 0.0582rs + 0.0079rs^2 - 0.0005rs^3

  G(q)::Float64 = A*(1 - exp(-B*(q/kF)^2))
end

function localfield_IU(nfe::NFE)::Function
  @attach(nfe, ρ)

  qF = fermiwavenumber(nfe)
  rs = electrondistance(nfe)
  α = (4/9π)^(1/3)

  b₀ = 0.0621814
  b₁ = 9.81379
  b₂ = 2.82224
  b₃ = 0.736411

  x = √rs
  rsΔEc = b₀ * (1 + b₁*x) / (1 + b₁*x + b₂*x^2 + b₃*x^3)
  rsΔEc′ = b₀ * ( b₁*(1 + b₁*x + b₂*x^2 + b₃*x^3) - (1 + b₁*x)*(b₁ + 2b₂*x + 3b₃*x^2) ) / (1 + b₁*x + b₂*x^2 + b₃*x^3)^2

  γ₀ = 1/4 - π*α/24 * rs^5 * (-3*rs^-4*rsΔEc + rs^-3 * -1/2/√rs * rsΔEc′)

  z = 4*(α*rs/π)^(1/2)
  I₁(z) = besseli(1, z)
  g₀ = 1/8*(z/I₁(z))^2

  A = 0.029
  B = 9/16*γ₀ - 3/64*(1 - g₀) - 16/15*A
  C = -3/4*γ₀ + 9/16*(1 - g₀) - 16/5*A

  function G(q::Real)::Float64
    Q = q / qF
    A*Q^4 + B*Q^2 + C + (A*Q^4 + (B + 8/3*A)*Q^2 - C) * (4 - Q^2)/4Q * log(abs((2+Q)/(2-Q)))
  end
end

#
# Wavenumber-energy characteristic
# Ref: W. A. Harrison: Elementary Electronic Structure Revised Edition (2004), 462
#
function wnechar(nfe::NFE)::Array{Function,2}
  ρ::Float64 = nfe.ρ
  z::Vector{Float64} = nfe.pseudo.z

  ω::Vector{Function} = formfactor(nfe)
  ϵ::Function = dielectric(nfe)
  G::Function = localfield(nfe)

  N::Int = length(z)
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    F(q)::Float64 = - q^2 / (8π*ρ) * ω[i](q) * ω[j](q) / (1 / (ϵ(q) - 1) + (1 - G(q)))

    ret[i,j] = ret[j,i] = F
  end

  ret
end

function coreradius(nfe::NFE)::Array{Float64,2}
  N::Int = ncomp(nfe)
  rc::Vector{Float64} = nfe.pseudo.rc

  [ (rc[i] + rc[j]) / 2 for i in 1:N, j in 1:N ]
end

function cutoffradius(nfe::NFE)::Array{Float64,2}
  10 * coreradius(nfe)
end

#
# Effective pair-potential between ions including full Coulomb and indirect parts
# Ref: W. A. Harrison: Elementary Electronic Structure Revised Edition (2004), 490
#
function pairpotential(nfe::NFE)::Array{Function,2}
  ρ::Float64 = nfe.ρ
  z::Vector{Float64} = nfe.pseudo.z

  F::Array{Function,2} = wnechar(nfe)

  N::Int = length(z)
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    function u(r)::Float64
      Ft(r)::Float64 = ∫(q -> F[i,j](q) * sin(q*r)/(q*r) * q^2, 0, Q_MAX)
      z[i]*z[j] / r + 1 / (π^2*ρ) * Ft(r)
    end

    ret[i,j] = ret[j,i] = u
  end

  ret
end

function pairpotential_minimizer(nfe::NFE)::Array{Float64,2}
  N::Int = length(nfe.z)
  rc::Vector{Float64} = nfe.pseudo.rc
  u::Array{Function,2} = pairpotential(nfe)

  ret = Array{Float64,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    r̄c = (rc[i] + rc[j]) / 2

    opt = Optim.optimize(u[i,j], 4*√r̄c, 8*√r̄c)
    ret[i,j] = ret[j,i] =  Optim.minimizer(opt)
  end

  ret
end

function pairpotential_derivative(nfe::NFE)::Array{Function,2}
  ρ::Float64 = nfe.ρ
  z::Vector{Float64} = nfe.pseudo.z

  F::Array{Function,2} = wnechar(nfe)

  N::Int = length(z)
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    function u′(r)::Float64
      Ft(r)::Float64 =
        ∫(q -> F[i,j](q) * (cos(q*r)/r - sin(q*r)/(q*r^2)) * q^2, 0, Q_MAX)
      - z[i]*z[j] / r^2 + 1 / (π^2*ρ) * Ft(r)
    end

    ret[i,j] = ret[j,i] = u′
  end

  ret
end

function entropy(nfe::NFE, ref::ReferenceSystem, T::Float64)::Float64
  z = nfe.z
  kF = fermiwavenumber(nfe)

  S = z*T*(π*kB/kF)^2 / kB
end

function internal(nfe::NFE, ref::ReferenceSystem)::Float64
  U_eg = internal_eg(nfe)
  U_es = internal_es(nfe, ref)
  U_pair = internal_pair(nfe, ref)

  U = U_eg + U_es + U_pair
end

# Internal energy of electron gas
function internal_eg(nfe::NFE)::Float64
  z::Float64 = nfe.z
  kF::Float64 = fermiwavenumber(nfe)

  E_PN = -0.0474 - 0.0155*log(kF) # Pines-Nozieres exchange-correlation

  U_e = z * (0.3kF^2 - 3/4π*kF + E_PN)
end

#
# Electrostatic energy between ions and electrons
# Ref: W. A. Harrison: Elementary Electronic Structure Revised Edition (2004), 490
#
function internal_es(nfe::NFE, ref::ReferenceSystem)::Float64
  N::Int = ncomp(ref)
  ρ::Float64 = nfe.ρ
  c::Vector{Float64} = composition(ref)

  F::Array{Function,2} = wnechar(nfe)
  U_ρ::Float64 = 0

  for i in 1:N
    U_ρ += 1 / (2π^2*ρ) * c[i] * ∫(q -> F[i,i](q)*q^2, 0, Q_MAX)
  end

  U_ρ
end

function internal_pair(nfe::NFE, ref::ReferenceSystem)::Float64
  internal(ref, nfe)
end
