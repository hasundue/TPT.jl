#
# Common functions for ReferenceSystem
#
function nndistance(ref::ReferenceSystem, g::Array{Function,2})::Array{Float64,2}
  N::Int = ncomp(ref)
  σ::Array{Float64,2} = hsdiameter(ref)
  d = Array{Float64,2}(N,N)

  for i in 1:N, j in 1:N
    i > j && continue

    opt = Optim.optimize(r -> -g[i,j](r), 0.5σ[i,j], 1.5σ[i,j])
    d[i,j] = d[j,i] = Optim.minimizer(opt)
  end

  return d
end

function nndistance(ref::ReferenceSystem)::Array{Float64,2}
  g::Array{Function,2} = paircorrelation(ref)
  return nndistance(ref, g)
end
