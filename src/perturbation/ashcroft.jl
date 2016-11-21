"""
ashcroft.jl

Ashcroft's empty-core pseudopotential
"""

immutable Ashcroft <: PseudoPotential
  z::Vector{Float64} # valence of core ion
  rc::Vector{Float64} # pseudopotential core radius
end

function Ashcroft(z::Real, rc::Real)
  Ashcroft([z], [rc])
end

function formfactor(nfe::NFE{Ashcroft})::Vector{Function}
  ρ = nfe.ρ

  N = length(nfe.pseudo.z)
  ω = Vector{Function}(N)

  for i = 1:N
    z = nfe.pseudo.z[i]
    rc = nfe.pseudo.rc[i]

    ω[i] = q -> -4π*ρ / q^2 * z * cos(q*rc)
  end

  return ω
end
