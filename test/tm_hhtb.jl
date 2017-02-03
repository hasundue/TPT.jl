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
r₀ = (p[:V] * 3/4π) .^ (1/3)
Nd = p[:Nd]
Ed = p[:Ed]
Wd = p[:Wd]

#
# Binary alloys
#
x = [0.5, 0.5]
for a in 1:N, b in 1:N
  a ≠ 6 && continue
  a == b && continue

  println("$(X[a])-$(X[b])")

  hhtb = TPT.HHTB(x, [ Nd[a], Nd[b] ], [ Ed[a], Ed[b] ], [ Wd[a], Wd[b] ], [ r₀[a], r₀[b] ])

  Emin, Emax = TPT.cutoffenergy(hhtb)
  Ef = TPT.fermienergy(hhtb)
  D = TPT.partialdensityofstate(hhtb)
  D_total = TPT.totaldensityofstate(hhtb)

  plot([ D[1], D[2], D_total ], Emin, Emax, labels=[ X[a] X[b] "total" ])
  vline!([Ef], label="E_F")
  xlabel!("E (eV)")
  ylabel!("DOS (states / eV atom)")
  png(joinpath(resdir, "$(a)-$(b)_$(X[a])-$(X[b])_D"))

  E_band = TPT.bandenergy(hhtb)
  E_site = TPT.onsiteenergy(hhtb)

  hhtb_a = TPT.HHTB(Nd[a], Ed[a], Wd[a], r₀[a])
  hhtb_b = TPT.HHTB(Nd[b], Ed[b], Wd[b], r₀[b])

  E_band_a = TPT.bandenergy(hhtb_a)
  E_site_a = TPT.onsiteenergy(hhtb_a)
  E_band_b = TPT.bandenergy(hhtb_b)
  E_site_b = TPT.onsiteenergy(hhtb_b)

  ΔE_band = E_band - x[1]*E_band_a - x[2]*E_band_b
  println("ΔE_band = $(ΔE_band)")

  ΔE_site = E_site - x[1]*E_site_a - x[2]*E_site_b
  println("ΔE_site = $(ΔE_site)")
end
