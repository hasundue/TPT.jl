"""
bretonnet_silbert.jl

Bretonnet-Silbert pseudopotential
"""

immutable BretonnetSilbert <: PseudoPotential
  rc::Vector{Float64} # pseudopotential core radius
  a::Vector{Float64} # softness parameters
end

function BretonnetSilbert(rc::Number, a::Number)
  BretonnetSilbert([rc], [a])
end

function formfactor(nfe::NFE{BretonnetSilbert})
  ρ = nfe.ρ

  N = length(nfe.z)
  ω = Array{Function}(N)

  for i = 1:N
    z = nfe.z[i]
    rc = nfe.pp.rc[i]
    a = nfe.pp.a[i]

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
