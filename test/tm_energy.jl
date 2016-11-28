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
res = Array{Vector{TPT.TPTSystem},2}(M,M)

Threads.@threads for b in 1:M
  # (a,b) = ind2sub((M,M), k)

  # a ≥ b && continue

  a = 5

  ΔFm[a,b] = Vector{Float64}(11)

  T = p[:T][a]
  σ = [ p[:σ][a], p[:σ][b] ]
  ρ = [ p[:ρ][a], p[:ρ][b] ]

  m = [ p[:m][a], p[:m][b] ]
  zs = [ p[:zs][a], p[:zs][b] ]
  rc = [ p[:rc][a], p[:rc][b] ]
  pa = [ p[:a][a], p[:a][b] ]
  zd = [ p[:zd][a], p[:zd][b] ]
  rd = [ p[:rd][a], p[:rd][b] ]

  F = Vector{Float64}(11)
  res[a,b] = Vector{TPT.TPTSystem}(11)

  # Pure components
  for i in 1:2
    ahs = TPT.AHS(ρ = ρ[i], σ = σ[i], approx = "RFA")
    wca = TPT.LWCA(ahs, T)

    pse = TPT.BretonnetSilbert(zs[i], rc[i], pa[i])
    nfe = TPT.NFE(ahs, pse)
    tb = TPT.WHTB(zd[i], rd[i])
    nfetb = TPT.NFETB(nfe, tb)

    sys = TPT.TPTSystem(wca, nfetb, m = m[i])

    j = i == 1 ? 1 : 11

    F[j] = TPT.helmholtz(sys)
    res[a,b][j] = sys
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

    res[a,b][i] = sys
    F[i] = TPT.helmholtz(sys)
  end

  for i in 1:11
    x₂ = (i-1) / 10
    ΔFm[a,b][i] = F[i] - (1-x₂)*F[1] - x₂*F[11]
  end
end

for a in 5, b in 1:M
  # a ≥ b && continue

  A = p[:X][a]
  B = p[:X][b]
  sys = res[a,b]::Vector{TPT.TPTSystem}

  default(xlabel = "x2")

  # Free energy of mixing
  plot(0:0.1:1.0, 2625.5 * ΔFm[a,b], ylabel = "dFmix (kJ/mol)")
  png(joinpath(resdir, "$(A)-$(B)_F"))

  # Effective hard-sphere diameter
  σ₁ = [ TPT.hsdiameter(sys[i].ref)[1,1] for i in 2:10 ]
  σ₁ = cat(1, [TPT.hsdiameter(sys[1].ref)[1,1]], σ₁, [NaN])
  σ₂ = [ TPT.hsdiameter(sys[i].ref)[2,2] for i in 2:10 ]
  σ₂ = cat(1, [NaN], σ₂, [TPT.hsdiameter(sys[2].ref)[2,2]])
  plot(0:0.1:1, [σ₁, σ₂], labels = ["1" "2"],
       ylabel = "Effective HS diameter (a.u.)")
  png(joinpath(resdir, "$(A)-$(B)_sigma"))
end

@testset "TM Energy" begin
end
