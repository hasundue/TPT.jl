"""
ashcroft.jl

Ashcroft's empty-core pseudopotential
"""

immutable Ashcroft <: PseudoPotential
  rc::Vector{Float64} # pseudopotential core radius
end

function Ashcroft(rc::Number)
  Ashcroft([rc])
end

function formfactor(nfe::NFE{Ashcroft})
  ρ = nfe.ρ

  N = length(nfe.z)
  ω = Array{Function}(N)

  for i = 1:N
    z = nfe.z[i]
    rc = nfe.pp.rc[i]

    ω[i] = q -> -4π*ρ / q^2 * z * cos(q*rc)
  end

  return ω
end
