"""
hausleitner.jl

Bond order tight-binding approach by Hausleitner and Hafner.

References:

* C. Hausleitner and J. Hafner: Phys. Rev. B, 45 (1992) 115-127.
"""

immutable HHTB <: TBPerturbation
  N::Int # number of components
  Z::Int # total coordination number on the Bethe lattice
  σ::Float64 # ordering parameter
  x::Vector{Float64} # composition
  p::Matrix{Float64} # probability
  Nd::Vector{Float64} # number of d-electrons
  Ed::Vector{Float64} # d-state energy
  Wd::Vector{Float64} # band width
  r₀::Vector{Float64} # atomic radius
end

function HHTB(Nd::Real, Ed::Real, Wd::Real, r₀::Real; Z::Int = 12)
  HHTB(1, Z, 0, ones(1), ones(1,1), [Nd], [Ed], [Wd], [r₀])
end

function HHTB(x::Vector{Float64}, Nd::Vector{Float64}, Ed::Vector{Float64}, Wd::Vector{Float64}, r₀::Vector{Float64}; Z::Int = 12, σ::Real = 0)
  N = length(x)
  p = [ x[j] for i in 1:N, j in 1:N ]
  HHTB(N, Z, σ, x, p, Nd, Ed, Wd, r₀)
end

ncomp(hhtb::HHTB) = hhtb.N

function coreradius(hhtb)::Matrix{Float64}
  @attach(hhtb, N, r₀)
  [ 0.5 * (r₀[i] + r₀[j]) / 2 for i in 1:N, j in 1:N ]
end

function cutoffradius(hhtb)::Matrix{Float64}
  @attach(hhtb, N, r₀)
  [ 10 * (r₀[i] + r₀[j]) / 2 for i in 1:N, j in 1:N ]
end

function cutoffenergy(hhtb)::Tuple{Float64,Float64}
  @attach(hhtb, Ed, Wd)
  Emin = minimum(Ed) - maximum(Wd)
  Emax = maximum(Ed) + maximum(Wd)
  (Emin, Emax)
end

function nntransferintegral(hhtb::HHTB)::Matrix{Float64}
  @attach(hhtb, Z, N, Wd)
  hd = Matrix{Float64}(N,N)
  for i in 1:N
    hd[i,i] = Wd[i] / 12
  end
  for i in 1:N, j in 1:N
    i == j && continue
    hd[i,j] = √( hd[i,i] * hd[j,j] )
  end
  hd
end

function transferintegral(hhtb::HHTB)::Matrix{Function}
  @attach(hhtb, N, Wd, r₀)
  ret = Matrix{Function}(N,N)
  η_ddσ = -6 * 2/5
  η_ddπ =  4 * 2/5
  η_ddδ = -1 * 2/5
  for i in 1:N
    h(r) = √( 1/5 * (η_ddσ^2 + 2η_ddπ^2 + 2η_ddδ^2) ) * Wd[i] * (r₀[i]/r)^5
    ret[i,i] = h
  end
  for i in 1:N, j in 1:N
    i == j && continue
    h(r) = √( ret[i,i](r) * ret[j,j](r) )
    ret[i,j] = h
  end
  ret
end

function transfermatrix(hhtb::HHTB)::Matrix{Function}
  @attach(hhtb, Z, N, x, Ed, p)
  h::Matrix{Float64} = nntransferintegral(hhtb)
  ret = Matrix{Function}(N,N)
  if N == 1
    function S(E::Real)
      A₂ = Z*h[1,1]
      A₁ = -(E - Ed[1])
      A₀ = h[1,1]
      roots = Polynomials.roots( Poly([A₀, A₁, A₂]) )
      roots[1]
    end
    ret[1,1] = S
    return ret
  end
  function S₁₁(E::Real)
    A₄ = (Z-1)^2 * ( p[1,1]^2*p[2,2] - p[1,1]*p[1,2]*p[2,1] )
    A₃ = (Z-1) * ( (p[1,2]*p[2,1] - 2p[1,1]*p[2,2]) * (E - Ed[1]) / h[1,1] + p[1,1]*p[1,2] * (E - Ed[2]) / h[2,2] )
    A₂ = p[2,2] * ( (E - Ed[1]) / h[1,1] )^2 - p[1,2] * (E - Ed[1])*(E - Ed[2]) / ( h[1,1] * h[2,2] ) + (Z-1) * (2*p[1,1]*p[2,2] + p[1,2]^2 - p[1,2]*p[2,1])
    A₁ = p[1,2]*(E - Ed[2])/h[2,2] - 2*p[2,2]*(E - Ed[1])/h[1,1]
    A₀ = p[2,2]
    roots::Vector{Number} = cardano(A₃, A₂, A₁, A₀)
    roots[2]
  end
  S₁₂(E) = h[1,2]/h[1,1]*S₁₁(E)
  S₂₁(E) = (-(Z-1)*p[1,1]*h[1,1]*S₁₁(E)^2 + (E-Ed[1])*S₁₁(E) - h[1,1]) / ((Z-1)*p[1,2]*h[1,2]*S₁₁(E))
  S₂₂(E) = h[2,2]/h[2,1]*S₂₁(E)
  [ S₁₁ S₁₂; S₂₁ S₂₂ ]
end

function selfenergy(hhtb::HHTB)::Vector{Function}
  @attach(hhtb, Z, N, x, p)
  h::Matrix{Float64} = nntransferintegral(hhtb)
  S::Matrix{Function} = transfermatrix(hhtb)
  [ Δ(E) = Z*sum( p[i,j]*h[i,j]*S[j,i](E) for j in 1:N ) for i in 1:N ]
end

function diagonalgreenfunction(hhtb::HHTB)::Vector{Function}
  @attach(hhtb, Z, N, Ed)
  Δ::Vector{Function} = selfenergy(hhtb)
  [ Gᵢᵢ(E) = 1 / (E - Ed[i] - Δ[i](E)) for i in 1:N ]
end

function offdiagonalgreenfunction(hhtb::HHTB)::Matrix{Function}
  @attach(hhtb, Z, N, Ed)
  h::Matrix{Float64} = nntransferintegral(hhtb)
  Δ::Vector{Function} = selfenergy(hhtb)
  S::Matrix{Function} = transfermatrix(hhtb)
  [ Gᵢⱼ(E) = h[i,j] / ( E - Ed[i] - (Z-1)/Z * Δ[i](E) ) / ( E - Ed[j] - Δ[j](E) ) for i in 1:N, j in 1:N ]
end

function partialdensityofstate(hhtb::HHTB)::Vector{Function}
  @attach(hhtb, N, Z, x, Ed)
  G::Vector{Function} = diagonalgreenfunction(hhtb)
  [ D(E) = 10*x[i]*1/π*imag(G[i](E)) for i in 1:N ]
end

function totaldensityofstate(hhtb::HHTB)::Function
  @attach(hhtb, N)
  Dᵢ = partialdensityofstate(hhtb)
  D(E) = sum(Dᵢ[i](E) for i in 1:N)
end

densityofstate(hhtb) = totaldensityofstate(hhtb)

function fermienergy(hhtb::HHTB)::Float64
  @attach(hhtb, N, x, Nd, Ed, Wd)
  D::Function = totaldensityofstate(hhtb)
  Emin, Emax = cutoffenergy(hhtb)
  fopt(E) = abs(sum(x .* Nd) - ∫(D, Emin, E))
  opt = Optim.optimize(fopt, Emin, Emax)
  Ef = Optim.minimizer(opt)
end

function bondorder(hhtb::HHTB)::Matrix{Float64}
  @attach(hhtb, N, Ed, Wd)
  G::Matrix{Function} = offdiagonalgreenfunction(hhtb)
  Ef::Float64 = fermienergy(hhtb)
  Emin, Emax = cutoffenergy(hhtb)
  Θ = Matrix{Float64}(N,N)
  for i in 1:N
    Θ[i,i] = 10 * 2/π * ∫(E -> imag(G[i,i](E)), Emin, Ef)
  end
  for i in 1:N, j in 1:N
    i == j && continue
    Θ[i,j] = 10 * 2/π * 1/2 * ∫(E -> imag(G[i,j](E) + G[j,i](E)), Emin, Ef)
  end
  Θ
end

function pairpotential_bond(hhtb::HHTB)::Matrix{Function}
  @attach(hhtb, N)
  h::Matrix{Function} = transferintegral(hhtb)
  Θ::Matrix{Float64} = bondorder(hhtb)
  [ u_bond(r) = h[i,j](r) * Θ[i,j] for i in 1:N, j in 1:N ]
end

function pairpotential_rep(hhtb::HHTB)::Matrix{Function}
  @attach(hhtb, N, Nd, Wd, r₀)
  [ u_rep(r) = 8/25*√(Nd[i]*Nd[j])*Wd[i]*Wd[j]*r₀[i]^5*r₀[j]^5 / r^8
    for i in 1:N, j in 1:N ]
end

function pairpotential(hhtb::HHTB)::Matrix{Function}
  @attach(hhtb, N)
  u_bond = pairpotential_bond(hhtb)
  u_rep = pairpotential_rep(hhtb)
  [ u(r) = u_rep[i,j](r) + u_bond[i,j](r) for i in 1:N, j in 1:N ]
end

function chargetransfer(hhtb::HHTB)::Vector{Float64}
  @attach(hhtb, N, x, Nd)
  Emin, Emax = cutoffenergy(hhtb)
  Ef::Float64 = fermienergy(hhtb)
  D::Vector{Function} = partialdensityofstate(hhtb)
  [ ∫(D[i], Emin, Ef) / x[i] - Nd[i] for i in 1:N ]
end

function bandenergy(hhtb::HHTB)::Float64
  D::Function = totaldensityofstate(hhtb)
  Ef::Float64 = fermienergy(hhtb)
  Emin, Emax = cutoffenergy(hhtb)
  E_band = ∫(E -> D(E)*E, Emin, Ef)
end

function bondenergy(hhtb::HHTB)::Float64
  @attach(hhtb, N, Ed)
  D::Vector{Function} = partialdensityofstate(hhtb)
  Ef::Float64 = fermienergy(hhtb)
  Emin, Emax = cutoffenergy(hhtb)
  E_bond = sum( ∫(E -> D[i](E)*(E - Ed[i]), Emin, Ef) for i in 1:N )
end

function onsiteenergy(hhtb::HHTB)::Float64
  @attach(hhtb, N, x, Nd, Ed)
  ΔNd = chargetransfer(hhtb)
  E_site = sum( x[i]*(Nd[i] + ΔNd[i])*Ed[i] for i in 1:N )
end

function internal(hhtb::HHTB, ref::ReferenceSystem)::Float64
  u_rep = pairpotential_rep(hhtb)
  U_band = bandenergy(hhtb)
  U_rep = internal(ref, hhtb, u_rep)
  U = U_band + U_rep
end
