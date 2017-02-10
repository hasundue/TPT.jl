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

# Convertion of parameters (Hartree units)
X = p[:X]
Vc = p[:Vc] / (5.292^3 * 1e-3) # a.u
r₀ = (Vc * 3/4π) .^ (1/3) # a.u.
Vl = p[:Vl] / (5.292^3 * 1e-3) # a.u.
ρ = 1 / Vl # a.u.
σ = p[:σ] / 5.291e-1 # Å -> a.u.
Ns = p[:Ns]
Nd = p[:Nd]
Ed = p[:Ed] / 27.21 # eV -> a.u.
Wd = p[:Wd] / 27.21 # eV -> a.u.
Rc = p[:Rc] / 5.292e-1 # Å -> a.u.

# some parameters
rmin, rmax = (2, 12)
umin, umax = (-0.04, 0.04)

#
# Puremetals
#
E_band = zeros(N)
E_bond = zeros(N)
E_rep = zeros(N)
E_nfe = zeros(N)

for i in 1:N
  ahs = TPT.AHS(σ = σ[i], ρ = ρ[i])
  pp = TPT.Ashcroft(Ns[i], Rc[i])
  nfe = TPT.NFE(ahs, pp, :IU)
  hhtb = TPT.HHTB(Nd[i], Ed[i], Wd[i], r₀[i])
  nfetb = TPT.NFETB(nfe, hhtb)

  u_nfe = TPT.pairpotential(nfe)[1,1]
  u_tb = TPT.pairpotential(hhtb)[1,1]
  u_tot = TPT.pairpotential(nfetb)[1,1]
  plot([u_nfe, u_tb, u_tot], 2, 12, ylims=(umin, umax), xl="r (a.u.)", yl="u (a.u.)", labels=["NFE" "TB" "total"])
  png(joinpath(resdir, "$(i)_$(X[i])_u"))

  E_band[i] = TPT.bandenergy(hhtb)
  E_bond[i] = TPT.internal(ahs, hhtb, TPT.pairpotential_bond(hhtb))
  E_rep[i] = TPT.internal(ahs, hhtb, TPT.pairpotential_rep(hhtb))
  E_nfe[i] = TPT.internal(ahs, hhtb, TPT.pairpotential(nfe))
end

plot(X, [E_bond, E_rep, E_nfe, E_bond + E_rep + E_nfe], labels=["bond" "rep" "nfe" "total"], m=:o, xl="X", yl="E (a.u. / atom)")
png(joinpath(resdir, "E"))

#
# Fe-based binary alloys
#
E_alloy = zeros(N)
ΔE_bond = zeros(N)
ΔE_rep = zeros(N)
ΔE_nfe = zeros(N)
x = [0.5, 0.5]
for i in 6, j in 1:N
  i == j && continue

  ρ_alloy = (ρ[i] + ρ[j]) / 2
  ahs = TPT.AHS(ρ = ρ_alloy, σ = [ σ[i], σ[j] ], c = [0.5, 0.5])
  pp = TPT.Ashcroft([ Ns[i], Ns[j] ], [ Rc[i], Rc[j] ])
  nfe = TPT.NFE(ahs, pp, :IU)
  hhtb = TPT.HHTB(x, [ Nd[i], Nd[j] ], [ Ed[i], Ed[j] ], [ Wd[i], Wd[j] ], [ r₀[i], r₀[j] ])
  nfetb = TPT.NFETB(nfe, hhtb)

  Emin, Emax = TPT.cutoffenergy(hhtb)
  Ef = TPT.fermienergy(hhtb)
  D = TPT.partialdensityofstate(hhtb)
  D_total = TPT.totaldensityofstate(hhtb)

  plot([ D[1], D[2], D_total ], Emin, Emax, labels=[ X[i] X[j] "total" ])
  vline!([Ef], label="E_F")
  xlabel!("E (a.u.)")
  ylabel!("DOS (states / a.u. atom)")
  png(joinpath(resdir, "$(i)-$(j)_$(X[i])-$(X[j])_D"))

  u = TPT.pairpotential(nfetb)
  plot([ u[1,1], u[1,2], u[2,2] ], rmin, rmax, labels=["1-1" "1-2" "2-2"], ylims=(umin, umax), xl="r (a.u.)", yl="u (a.u.)")
  png(joinpath(resdir, "$(i)-$(j)_$(X[i])-$(X[j])_u"))

  E_rep_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential_rep(hhtb))
  E_bond_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential_bond(hhtb))
  E_nfe_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential(nfe))
  E_alloy[j] = E_rep_alloy + E_bond_alloy + E_nfe_alloy

  ΔE_bond[j] = E_bond_alloy - x[1]*E_bond[i] - x[2]*E_bond[j]
  ΔE_rep[j] = E_rep_alloy - x[1]*E_rep[i] - x[2]*E_rep[j]
  ΔE_nfe[j] = E_nfe_alloy - x[1]*E_nfe[i] - x[2]*E_nfe[j]
end

plot(X, E_alloy, label="total", m=:o, xl="Fe-X", yl="E (a.u. / atom)")
png(joinpath(resdir, "E_alloy"))

plot(X, [ΔE_bond, ΔE_rep, ΔE_nfe, ΔE_bond + ΔE_rep + ΔE_nfe], labels=["bond" "rep" "nfe" "total"], m=:o, xl="Fe-X", yl="Emix (a.u. / atom)")
png(joinpath(resdir, "Emix"))

@testset "TM HHTB" begin
end
