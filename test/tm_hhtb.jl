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
V = p[:V] / (5.292^3 * 1e-3) # a.u.
r₀ = (V * 3/4π) .^ (1/3) # a.u.
ρ = 1 / V * 0.9 # a.u.
σ = p[:σ] / 5.291 * 10 # a.u.
Ns = p[:Ns]
Nd = p[:Nd]
Ed = p[:Ed] / 27.21 * 2 # Ry
Wd = p[:Wd] / 27.21 * 2 # Ry
Rc = p[:Rc] / 5.292e-1 # a.u.

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
  nfe = TPT.NFE(ahs, pp)
  hhtb = TPT.HHTB(Nd[i], Ed[i], Wd[i], r₀[i])
  nfetb = TPT.NFETB(nfe, hhtb)
  # lwca = TPT.LWCA(ahs, 1873, struct=:full)
  # sys = TPT.TPTSystem(lwca, nfetb)

  # u_rep = TPT.pairpotential_rep(hhtb)
  u_nfe = TPT.pairpotential(nfe)[1,1]
  u_tb = TPT.pairpotential(hhtb)[1,1]
  u_tot = TPT.pairpotential(nfetb)[1,1]
  plot(u_nfe, 2, 12, ylims=(-0.2, 0.3), xl="r (a.u.)", yl="u (Ry)", label="NFE")
  plot!(u_tb, 2, 12, ylims=(-0.2, 0.3), xl="r (a.u.)", yl="u (Ry)", label="TB")
  plot!(u_tot, 2, 12, ylims=(-0.2, 0.3), xl="r (a.u.)", yl="u (Ry)", label="total")
  png(joinpath(resdir, "$(i)_$(X[i])_u"))

  g_ahs = TPT.paircorrelation(ahs)[1,1]
  g_wca = TPT.paircorrelation(sys)[1,1]
  σ_wca = TPT.hsdiameter(sys.ref)
  plot([g_ahs, g_wca], 2, 12, xl="r (a.u.)", yl="g(r)", labels=["AHS" "WCA"])
  vline!([σ_wca], label="eff.HSD")
  png(joinpath(resdir, "$(i)_$(X[i])_g"))

  E_band[i] = TPT.bandenergy(hhtb)
  # E_bond[i] = TPT.bondenergy(hhtb)
  E_bond[i] = TPT.internal(ahs, hhtb, TPT.pairpotential_bond(hhtb))
  E_rep[i] = TPT.internal(ahs, hhtb, TPT.pairpotential_rep(hhtb))
  E_nfe[i] = TPT.internal(ahs, hhtb, TPT.pairpotential(nfe)) * 2
end

plot(X, [E_bond, E_rep, E_nfe, E_bond + E_rep + E_nfe], labels=["bond" "rep" "nfe" "total"], m=:o, xl="X", yl="E (Ry / atom)")
png(joinpath(resdir, "E"))

#
# Binary alloys
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
  nfe = TPT.NFE(ahs, pp)
  hhtb = TPT.HHTB(x, [ Nd[i], Nd[j] ], [ Ed[i], Ed[j] ], [ Wd[i], Wd[j] ], [ r₀[i], r₀[j] ])
  nfetb = TPT.NFETB(nfe, hhtb)

  Emin, Emax = TPT.cutoffenergy(hhtb)
  Ef = TPT.fermienergy(hhtb)
  D = TPT.partialdensityofstate(hhtb)
  D_total = TPT.totaldensityofstate(hhtb)

  plot([ D[1], D[2], D_total ], Emin, Emax, labels=[ X[i] X[j] "total" ])
  vline!([Ef], label="E_F")
  xlabel!("E (Ry)")
  ylabel!("DOS (states / Ry atom)")
  png(joinpath(resdir, "$(i)-$(j)_$(X[i])-$(X[j])_D"))

  u = TPT.pairpotential(nfetb)
  plot([ u[1,1], u[1,2], u[2,2] ], 2, 8, labels=["1-1" "1-2" "2-2"], ylims=(-0.2, 0.3), xl="r (a.u.)", yl="u (Ry)")
  png(joinpath(resdir, "$(i)-$(j)_$(X[i])-$(X[j])_u"))

  E_rep_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential_rep(hhtb))
  # E_band_alloy = TPT.bandenergy(hhtb)
  E_bond_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential_bond(hhtb))
  E_nfe_alloy = TPT.internal(ahs, hhtb, TPT.pairpotential(nfe)) * 2
  E_alloy[j] = E_rep_alloy + E_bond_alloy + E_nfe_alloy

  # ΔE_bond[j] = E_bond_alloy - x[1]*E_bond[i] - x[2]*E_bond[j]
  ΔE_bond[j] = E_bond_alloy - x[1]*E_bond[i] - x[2]*E_bond[j]
  ΔE_rep[j] = E_rep_alloy - x[1]*E_rep[i] - x[2]*E_rep[j]
  ΔE_nfe[j] = E_nfe_alloy - x[1]*E_nfe[i] - x[2]*E_nfe[j]
end

plot(X, E_alloy, label="total", m=:o, xl="Fe-X", yl="E (Ry / atom)")
png(joinpath(resdir, "E_alloy"))

plot(X, [ΔE_bond, ΔE_rep, ΔE_nfe, ΔE_bond + ΔE_rep + ΔE_nfe], labels=["bond" "rep" "nfe" "total"], m=:o, xl="Fe-X", yl="Emix (Ry / atom)")
png(joinpath(resdir, "Emix"))

@testset "TM HHTB" begin
end
