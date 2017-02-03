import TPT
using Base.Test

using DataFrames
using Plots; pyplot()

println("--- TM HHTB ---")

# Load elemental parameters
p = readtable(joinpath("data", "parameters", "tm_hhtb.csv"))
N = size(p, 1) # number of elements

#
# Prepare a directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_hhtb")
!isdir(resdir) && mkdir(resdir)

# Convertion of parameters
X = p[:X]
r₀ = (p[:V] * 3/4π) .^ (1/3) # Å
ρ = 1/p[:V] # Å
σ = p[:σ]
Nd = p[:Nd]
Ed = p[:Ed] # eV
Wd = p[:Wd] # eV

#
# Puremetals
#
E_band = zeros(N)
E_rep = zeros(N)
for i in 1:N
  ahs = TPT.AHS(σ = σ[i], ρ = ρ[i])
  hhtb = TPT.HHTB(Nd[i], Ed[i], Wd[i], r₀[i])
  # u_rep = TPT.pairpotential_rep(hhtb)
  u = TPT.pairpotential(hhtb)[1,1]
  plot(u, 2, 8, ylims=(-0.3, 0.1))
  png(joinpath(resdir, "$(i)_$(X[i])_u"))
  # E_band[i] = TPT.bandenergy(hhtb)
  E_band[i] = TPT.internal(ahs, hhtb, TPT.pairpotential_bond(hhtb))
  E_rep[i] = TPT.internal(ahs, hhtb, TPT.pairpotential_rep(hhtb))
end

plot(X, [E_band, E_rep, E_band + E_rep], labels=["band" "rep" "total"], m=:o, xl="X", yl="E (eV / atom)")
png(joinpath(resdir, "E"))

#
# Binary alloys
#
ΔE_band = zeros(N)
ΔE_rep = zeros(N)
x = [0.5, 0.5]
for i in 6, j in 1:N
  i == j && continue

  hhtb = TPT.HHTB(x, [ Nd[i], Nd[j] ], [ Ed[i], Ed[j] ], [ Wd[i], Wd[j] ], [ r₀[i], r₀[j] ])

  Emin, Emax = TPT.cutoffenergy(hhtb)
  Ef = TPT.fermienergy(hhtb)
  D = TPT.partialdensityofstate(hhtb)
  D_total = TPT.totaldensityofstate(hhtb)

  plot([ D[1], D[2], D_total ], Emin, Emax, labels=[ X[i] X[j] "total" ])
  vline!([Ef], label="E_F")
  xlabel!("E (eV)")
  ylabel!("DOS (states / eV atom)")
  png(joinpath(resdir, "$(i)-$(j)_$(X[i])-$(X[j])_D"))

  u = TPT.pairpotential(hhtb)
  plot([ u[1,1], u[1,2], u[2,2] ], 2, 8, labels=["1-1" "1-2" "2-2"], ylims=(-0.3, 0.1))
  png(joinpath(resdir, "$(i)-$(j)_$(X[i])-$(X[j])_u"))

  ρ_alloy = (ρ[i] + ρ[j]) / 2
  ahs = TPT.AHS(ρ = ρ_alloy, σ = [ σ[i], σ[j] ], c = [0.5, 0.5])
  E_rep_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential_rep(hhtb))

  # E_band_alloy = TPT.bandenergy(hhtb)
  E_band_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential_bond(hhtb))

  ΔE_band[j] = E_band_alloy - x[1]*E_band[i] - x[2]*E_band[j]
  ΔE_rep[j] = E_rep_alloy - x[1]*E_rep[i] - x[2]*E_rep[j]
end

plot(X, [ΔE_band, ΔE_rep, ΔE_band + ΔE_rep], labels=["band" "rep" "total"], m=:o, xl="Fe-X", yl="Emix (eV / atom)")
# plot(X, ΔE_rep, label=:none, m=:o)
png(joinpath(resdir, "Emix"))
