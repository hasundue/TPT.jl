"""
bretonnet_silbert.jl

Bretonnet-Silbert pseudopotential
"""

immutable BretonnetSilbert <: PseudoPotential
  z::Vector{Float64} # valence of core ion
  rc::Vector{Float64} # pseudopotential core radius
  a::Vector{Float64} # softness parameters
end

function BretonnetSilbert(z::Real, rc::Real, a::Real)
  BretonnetSilbert([z], [rc], [a])
end

function formfactor(nfe::NFE{BretonnetSilbert})::Vector{Function}
  ρ::Float64 = nfe.ρ

  N::Int = length(nfe.pseudo.z)
  ret = Vector{Function}(N)

  for i in 1:N
    z::Float64 = nfe.pseudo.z[i]
    rc::Float64 = nfe.pseudo.rc[i]
    a::Float64 = nfe.pseudo.a[i]

    B₁::Float64 = (z/rc) * (1 - 2a/rc) * exp(rc/a)
    B₂::Float64 = (2z/rc) * (a/rc - 1) * exp(0.5rc/a)

    function ω(q)::Float64
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

    ret[i] = ω
  end

  return ret
end
