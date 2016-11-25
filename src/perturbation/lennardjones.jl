immutable LennardJones <: Perturbation
  ϵ::Array{Float64,1}
  σ::Array{Float64,1}
end

LennardJones(ϵ::Real, σ::Real) = LennardJones([ϵ], [σ])

ncomp(lj::LennardJones)::Int = length(lj.ϵ)

function pairpotential(lj::LennardJones)::Array{Function,2}
  @attach(lj, ϵ, σ)

  N::Int = ncomp(lj)
  ret = Array{Function,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    ϵᵢⱼ = √(ϵ[i]*ϵ[j])
    σᵢⱼ = (σ[i] + σ[j]) / 2

    function u(r)::Float64
      # if r < 2^(1/6) * σᵢⱼ
        4ϵᵢⱼ * ((σᵢⱼ/r)^12 - (σᵢⱼ/r)^6)
      # else
        # 0
      # end
    end

    ret[i,j] = ret[j,i] = u
  end

  return ret
end
