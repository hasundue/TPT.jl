"""
abstractwca.jl

Common datatypes and functions for WCA and LWCA

References:
"""

abstract AbstractWCA <: DependentReferenceSystem
abstract AbstractOptimizedWCA <: AbstractWCA

ncomp(wca::AbstractWCA)::Int = ncomp(wca.trial)
numberdensity(wca::AbstractWCA)::Vector{Float64} = numberdensity(wca.trial)
totalnumberdensity(wca::AbstractWCA)::Float64 = totalnumberdensity(wca.trial)
composition(wca::AbstractWCA)::Vector{Float64} = composition(wca.trial)
temperature(wca::AbstractWCA)::Float64 = wca.temp
hsdiameter(wca::AbstractWCA)::Array{Float64,2} = hsdiameter(wca.trial)

function emptyradius(wca::AbstractOptimizedWCA)::Array{Float64,2}
  N::Int = ncomp(wca)
  σ::Array{Float64,2} = hsdiameter(wca)
  B::Array{Function,2} = blipfunction(wca)

  ret = Array{Float64,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    opt = NLopt.Opt(:LN_BOBYQA, 1)
    NLopt.min_objective!(opt, (r, g) -> abs(B[i,j](r[1])))
    NLopt.lower_bounds!(opt, [0.])
    NLopt.upper_bounds!(opt, [σ[i,j]])
    NLopt.stopval!(opt, eps(Float64))

    (fmin, r₀, res) = NLopt.optimize(opt, [σ[i,j]])

    ret[i,j] = ret[j,i] = r₀[1]
  end

  return ret
end
