import TPT

using Base.Test

using DataFrames
using Optim
using Dierckx
using Plots; pyplot()

println("--- TM Energy ---")

#
# setting up the directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_energy")
!isdir(resdir) && mkdir(resdir)

# Elemental parameters
p = readtable(joinpath("data", "parameters", "tm_optim.csv"))

# number of elements available
M = size(p, 1)

# Results
ΔFm = Array{Vector{Float64},2}(M,M)

Threads.@threads for b in 1:M
  # (a,b) = ind2sub((M,M), k)

  # a ≥ b && continue

  a = 5
  a == b && continue

  ΔFm[a,b] = Vector{Float64}(11)

  T = (p[:T][a] + p[:T][b]) / 2
  σ = [ p[:σ][a], p[:σ][b] ]
  ρ = [ p[:ρ][a], p[:ρ][b] ]

  m = [ p[:m][a], p[:m][b] ]
  zs = [ p[:zs][a], p[:zs][b] ]
  rc = [ p[:rc][a], p[:rc][b] ]
  pa = [ p[:a][a], p[:a][b] ]
  zd = [ p[:zd][a], p[:zd][b] ]
  rd = [ p[:rd][a], p[:rd][b] ]

  F = Vector{Float64}(11)

  # Pure components
  for i in 1:2
    ahs = TPT.AHS(ρ = ρ[i], σ = σ[i], approx = "RFA")
    wca = TPT.LWCA(ahs, T)

    pse = TPT.BretonnetSilbert(zs[i], rc[i], pa[i])
    nfe = TPT.NFE(ahs, pse)
    tb = TPT.WHTB(zd, rd)
    nfetb = TPT.NFETB(nfe, tb)

    sys = TPT.TPTSystem(wca, nfetb, m = m[i])

    if i == 1
      F[1] = TPT.helmholtz(sys)
    else
      F[11] = TPT.helmholtz(sys)
    end
  end

  # Alloys
  for i in 2:10
    x₂ = (i-1) / 10

    ρ₀ = (1-x₂) * ρ[1] + x₂ * ρ[2]

    ahs = TPT.AHS(ρ = ρ₀, σ = σ, c = [1-x₂, x₂], approx = "RFA")
    wca = TPT.LWCA(ahs, T)

    pse = TPT.BretonnetSilbert(zs, rc, pa)
    nfe = TPT.NFE(ahs, pse)
    tb = TPT.WHTB(zd, rd)
    nfetb = TPT.NFETB(nfe, tb)

    sys = TPT.TPTSystem(wca, nfetb, m = m)

    F[i] = TPT.helmholtz(sys)
  end

  for i in 1:11
    x₂ = (i-1) / 10
    ΔFm[a,b][i] = F[i] - (1-x₂)*F[1] - x₂*F[11]
  end
end

default(xlabel = "x2", ylabel = "dFmix (kJ/mol)")

for a in 1:M, b in 1:M
  a ≥ b && continue

  A = p[:X][a]
  B = p[:X][b]

  plot(0:0.1:1.0, 2625.5 * ΔFm[a,b])
  png(joinpath(resdir, "$(A)-$(B)"))
end

@testset "TM Energy" begin
end
