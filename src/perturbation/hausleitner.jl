"""
hausleitner.jl

Bond order tight-binding approach by Hausleitner and Hafner.

References:

* C. Hausleitner and J. Hafner: Phys. Rev. B, 45 (1992) 115-127.
"""

immutable BOTB <: TBPerturbation
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

function BOTB(Nd::Real, Ed::Real, Wd::Real, r₀::Real; Z::Int = 12)
  BOTB(1, Z, 0, ones(1), ones(1,1), [Nd], [Ed], [Wd], [r₀])
end

function BOTB(x::Vector{Float64}, Nd::Vector{Float64}, Ed::Vector{Float64}, Wd::Vector{Float64}, r₀::Vector{Float64}; Z::Int = 12, σ::Real = 0)
  N = length(x)
  p = [ x[j] for i in 1:N, j in 1:N ]
  BOTB(N, Z, σ, x, p, Nd, Ed, Wd, r₀)
end

ncomp(botb::BOTB) = botb.N

function nntransferintegral(botb::BOTB)::Matrix{Float64}
  @attach(botb, Z, N, Wd)
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

function transferintegral(botb::BOTB)::Matrix{Function}
  @attach(botb, N, Wd, r₀)
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

function transfermatrix1(botb::BOTB)::Function
  @attach(botb, Z, Ed)
  h::Matrix{Float64} = nntransferintegral(botb)
  function S(E::Real)
    A₂ = Z*h[1,1]
    A₁ = -(E - Ed[1])
    A₀ = h[1,1]
    roots = Polynomials.roots( Poly([A₀, A₁, A₂]) )
    roots[1]
  end
end

function selfenergy1(botb::BOTB)::Function
  @attach(botb, Z)
  h::Matrix{Float64} = nntransferintegral(botb)
  S::Function = transfermatrix1(botb)

  Δ(E) = Z*h[1,1]*S(E)
end

function greenfunction1(botb::BOTB)::Function
  @attach(botb, Z, Ed)
  Δ = selfenergy1(botb)
  G(E) = 1 / (E - Ed[1] - Δ(E))
end

function transfermatrix(botb::BOTB)
  @attach(botb, Z, N, x, Ed, p)
  if N == 1
    return transfermatrix1(botb::BOTB)
  end
  h::Matrix{Float64} = nntransferintegral(botb)
  ret = Matrix{Function}(N,N)
  function S₁₁(E::Real)
    A₄ = (Z-1)^2 * ( p[1,1]^2*p[2,2] - p[1,1]*p[1,2]*p[2,1] )
    A₃ = (Z-1) * ( (p[1,2]*p[2,1] - 2p[1,1]*p[2,2]) * (E - Ed[1]) / h[1,1] + p[1,1]*p[1,2] * (E - Ed[2]) / h[2,2] )
    A₂ = p[2,2] * ( (E - Ed[1]) / h[1,1] )^2 - p[1,2] * (E - Ed[1])*(E - Ed[2]) / ( h[1,1] * h[2,2] ) + (Z-1) * (2*p[1,1]*p[2,2] + p[1,2]^2 - p[1,2]*p[2,1])
    A₁ = p[1,2]*(E - Ed[2])/h[2,2] - 2*p[2,2]*(E - Ed[1])/h[1,1]
    A₀ = p[2,2]
    S::Vector{Number} = cardano(A₃, A₂, A₁, A₀)
    S[2]
  end
  S₁₂(E) = h[1,2]/h[1,1]*S₁₁(E)
  S₂₁(E) = (-(Z-1)*p[1,1]*h[1,1]*S₁₁(E)^2 + (E-Ed[1])*S₁₁(E) - h[1,1]) / ((Z-1)*p[1,2]*h[1,2]*S₁₁(E))
  S₂₂(E) = h[2,2]/h[2,1]*S₂₁(E)
  [ S₁₁ S₁₂; S₂₁ S₂₂ ]
end

function selfenergy(botb::BOTB)::Vector{Function}
  @attach(botb, Z, N, x, p)
  h::Matrix{Float64} = nntransferintegral(botb)
  S::Matrix{Function} = transfermatrix(botb)
  [ Δ₁(E) = Z * ( p[1,1]*h[1,1]*S[1,1](E) + p[1,2]*h[1,2]*S[2,1](E) ),
    Δ₂(E) = Z * ( p[2,2]*h[2,2]*S[2,2](E) + p[2,1]*h[2,1]*S[1,2](E) ) ]
end

function diagonalgreenfunction(botb::BOTB)::Vector{Function}
  @attach(botb, Z, N, Ed)
  Δ::Vector{Function} = selfenergy(botb)
  [ Gᵢᵢ(E) = 1 / (E - Ed[i] - Δ[i](E)) for i in 1:N ]
end

function offdiagonalgreenfunction(botb::BOTB)::Matrix{Function}
  @attach(botb, Z, N, Ed)
  h::Matrix{Float64} = nntransferintegral(botb)
  Δ::Vector{Function} = selfenergy(botb)
  S::Matrix{Function} = transfermatrix(botb)
  [ Gᵢⱼ(E) = h[i,j] / ( E - Ed[i] - (Z-1)/Z * Δ[i](E) ) / ( E - Ed[j] - Δ[j](E) ) for i in 1:N, j in 1:N ]
end

function partialdensityofstate(botb::BOTB)::Vector{Function}
  @attach(botb, N, Z, x, Ed)
  G::Vector{Function} = diagonalgreenfunction(botb)
  [ D(E) = 10*x[i]*1/π*imag(G[i](E)) for i in 1:N ]
end

function totaldensityofstate(botb::BOTB)::Function
  @attach(botb, N)
  Dᵢ = partialdensityofstate(botb)
  D(E) = sum(Dᵢ[i](E) for i in 1:N)
end

densityofstate(botb) = totaldensityofstate(botb)

function fermienergy(botb::BOTB)::Float64
  @attach(botb, N, x, Nd, Ed, Wd)
  D::Function = totaldensityofstate(botb)
  Emin = minimum(Ed) - maximum(Wd)
  Emax = maximum(Ed) + maximum(Wd)
  fopt(E) = abs(sum(x .* Nd) - ∫(D, Emin, E))
  opt = Optim.optimize(fopt, Emin, Emax)
  Ef = Optim.minimizer(opt)
end

function bondorder(botb::BOTB)::Matrix{Float64}
  @attach(botb, N, Ed, Wd)
  G::Matrix{Function} = offdiagonalgreenfunction(botb)
  Ef::Float64 = fermienergy(botb::BOTB)
  Emin = minimum(Ed) - maximum(Wd)
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

function pairpotential_bond(botb::BOTB)::Matrix{Function}
  @attach(botb, N)
  h::Matrix{Function} = transferintegral(botb)
  Θ::Matrix{Float64} = bondorder(botb)
  [ u_bond(r) = h[i,j](r) * Θ[i,j] for i in 1:N, j in 1:N ]
end

function pairpotential_rep(botb::BOTB)::Matrix{Function}
  @attach(botb, N, Nd, Wd, r₀)
  [ u_rep(r) = 8/25*√(Nd[i]*Nd[j])*Wd[i]*Wd[j]*r₀[i]^5*r₀[j]^5 / r^8
    for i in 1:N, j in 1:N ]
end

function pairpotential(botb::BOTB)::Matrix{Function}
  @attach(botb, N)
  u_bond = pairpotential_bond(botb)
  u_rep = pairpotential_rep(botb)
  [ u(r) = u_rep[i,j](r) + u_bond[i,j](r) for i in 1:N, j in 1:N ]
end
