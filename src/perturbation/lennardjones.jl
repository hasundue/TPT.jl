immutable LennardJones <: Perturbation
  ϵ::Array{Float64,1}
  σ::Array{Float64,1}
end

LennardJones(ϵ::Real, σ::Real) = LennardJones([ϵ], [σ])

ncomp(lj::LennardJones)::Int = length(lj.ϵ)

function coreradius(lj::LennardJones)::Array{Float64,2}
  0.5 * pairpotential_minimizer(lj)
end

function cutoffradius(lj::LennardJones)::Array{Float64,2}
  5 * pairpotential_minimizer(lj)
end

function hsdiameter_estimate(lj::LennardJones)::Array{Float64,2}
  0.9 * pairpotential_minimizer(lj)
end

function pairpotential(lj::LennardJones)::Array{Function,2}
  @attach(lj, ϵ, σ)

  N::Int = ncomp(lj)
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    ϵᵢⱼ = √(ϵ[i]*ϵ[j])
    σᵢⱼ = (σ[i] + σ[j]) / 2

    u(r::Real)::Float64 = 4ϵᵢⱼ * ((σᵢⱼ/r)^12 - (σᵢⱼ/r)^6)

    ret[i,j] = ret[j,i] = u
  end

  return ret
end

function pairpotential_minimizer(lj::LennardJones)::Array{Float64,2}
  @attach(lj, ϵ, σ)

  N::Int = ncomp(lj)

  [ 2^(1/6) * (σ[i] + σ[j]) / 2 for i in 1:N, j in 1:N ]
end

function pairpotential_derivative(lj::LennardJones)::Array{Function,2}
  @attach(lj, ϵ, σ)

  N::Int = ncomp(lj)
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    ϵᵢⱼ = √(ϵ[i]*ϵ[j])
    σᵢⱼ = (σ[i] + σ[j]) / 2

    u′(r::Real)::Float64 = 4ϵᵢⱼ * (-12σᵢⱼ^12 / r^13 + 6σᵢⱼ^6 / r^7)

    ret[i,j] = ret[j,i] = u′
  end

  return ret
end
