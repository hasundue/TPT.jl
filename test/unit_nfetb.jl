import TPT
using Base.Test
using Plots; pyplot()

println("--- Unit NFETB ---")

#
# Liquid Fe
# The parameters are taken from: Bretonnet et al.: J. Phys.: Condens.
# Matter, 7 (1995), 3803-3815
# a is modified so that WCA can produce realistic g(r).
#
Ω = 89.35698
ρ = 1 / Ω
T = 1833.
zs = 1.4
rc = 1.540
a = 0.333
zd = 6.6
rd = 1.512

ahs = TPT.AHS(σ = 4.25, ρ = ρ, approx=:RFA)
lwca = TPT.LWCA(ahs, T)
wca = TPT.WCA(ahs, T)

nfe = TPT.NFE(ρ, zs, TPT.BretonnetSilbert(zs, rc, a))
tb = TPT.WHTB(zd, rd, version=:original)
nfetb = TPT.NFETB(nfe, tb)

sys1 = TPT.TPTSystem(ahs, nfetb, T = 1833., m = 1.018e+5)
sys2 = TPT.TPTSystem(wca, nfetb, m = 1.018e+5)
sys3 = TPT.TPTSystem(lwca, nfetb, m = 1.018e+5)


σ_estim = TPT.hsdiameter_estimate(nfetb, 1883.)[1,1]
σ1 = TPT.hsdiameter(sys1.ref)[1,1]
σ2 = TPT.hsdiameter(sys2.ref)[1,1]
σ3 = TPT.hsdiameter(sys3.ref)[1,1]

u_nfe = TPT.pairpotential(nfe)[1,1]
u_tb = TPT.pairpotential(tb)[1,1]
u_tot = TPT.pairpotential(nfetb)[1,1]
u′_tot = TPT.pairpotential_derivative(nfetb)[1,1]
rmin = TPT.pairpotential_minimizer(nfetb)[1,1]

plot([u_nfe, u_tb, u_tot, u′_tot], 2, 10, ylims = (-0.06, 0.1),
     labels = ["NFE" "TB" "Total" "derivative"])

g1 = TPT.paircorrelation(sys1)[1,1]
g2 = TPT.paircorrelation(sys2)[1,1]
g3 = TPT.paircorrelation(sys3)[1,1]
plot([g1, g2, g3], 3, 7, labels = ["AHS" "WCA" "LWCA"], ylims=:auto)

B2 = TPT.blipfunction(sys2.ref)[1,1]
B3 = TPT.blipfunction(sys3.ref)[1,1]
plot([r -> B2(r), r -> B3(r)], 3.0, 5.5, ylims=:auto, labels=["WCA" "LWCA"])

S1 = TPT.structurefactor(sys1.ref)[1,1]
S2 = TPT.structurefactor(sys2.ref)[1,1]
S3 = TPT.structurefactor(sys3.ref)[1,1]
plot([S1, S2, S3], 0.1, 6, labels = ["AHS" "WCA" "LWCA"], ylims=:auto)


# Internal energy of reference systems
U1_ref = TPT.internal(sys1.ref)
U2_ref = TPT.internal(sys2.ref)

# Internal energy from pairwise-contribution of NFE
U1_pair = TPT.internal_pair(sys1.pert.nfe, sys1.ref)
U2_pair = TPT.internal_pair(sys2.pert.nfe, sys2.ref)

# Internal energy from electrostatic contribution of NFE
U1_es = TPT.internal_es(sys1.pert.nfe, sys1.ref)
U2_es = TPT.internal_es(sys2.pert.nfe, sys2.ref)

# Internal energy of electron gas
U1_eg = TPT.internal_eg(sys1.pert.nfe)
U2_eg = TPT.internal_eg(sys2.pert.nfe)

# Total internal energy of NFE
U1_nfe = TPT.internal(sys1.pert.nfe, sys1.ref)
U2_nfe = TPT.internal(sys2.pert.nfe, sys2.ref)

# Band width
Wd1 = TPT.bandwidth(tb, sys1.ref)
Wd2 = TPT.bandwidth(tb, sys2.ref)

# Bonding energy of TB
U1_bond = TPT.internal_band(sys1.pert.tb, sys1.ref, Wd1)
U2_bond = TPT.internal_band(sys2.pert.tb, sys2.ref, Wd2)

# Pair-wise and repusive contribution
U1_rep = TPT.internal_rep(sys1.pert.tb, sys1.ref)
U2_rep = TPT.internal_rep(sys2.pert.tb, sys2.ref)

# Total TB energy
U1_tb = TPT.internal(sys1.pert.tb, sys1.ref)
U2_tb = TPT.internal(sys2.pert.tb, sys2.ref)

# Total internal energy
U1 = TPT.internal(sys1)
U2 = TPT.internal(sys2)

# Total entropy
S1 = TPT.entropy(sys1)
S2 = TPT.entropy(sys2)

# Total free energy
F1 = TPT.helmholtz(sys1)
F2 = TPT.helmholtz(sys2)

@testset "Unit NFETB" begin
  @test isapprox(U1_ref, 0.00, atol=1e-2)
  @test isapprox(U2_ref, 0.00861, atol=1e-5)

  @test isapprox(U1_pair, 0.0238, atol=1e-4)
  @test isapprox(U2_pair, 0.0278, atol=1e-4)
  @test isapprox(U1_es, -0.640, atol=1e-3)
  @test isapprox(U2_es, -0.640, atol=1e-3)
  @test isapprox(U1_nfe, -0.684, atol=1e-3)
  @test isapprox(U2_nfe, -0.680, atol=1e-3)

  @test isapprox(Wd1, 0.155, atol=1e-3)
  @test isapprox(Wd2, 0.170, atol=1e-3)
  @test isapprox(U1_bond, -0.174, atol=1e-3)
  @test isapprox(U2_bond, -0.191, atol=1e-3)
  @test isapprox(U1_rep, 0.0394, atol=1e-4)
  @test isapprox(U2_rep, 0.0420, atol=1e-4)
  @test isapprox(U1_tb, -0.135, atol=1e-3)
  @test isapprox(U2_tb, -0.149, atol=1e-3)

  @test isapprox(U1, -0.819, atol=1e-3)
  @test isapprox(U2, -0.829, atol=1e-3)
  @test isapprox(S1, 11.2, atol=1e-1)
  @test isapprox(S2, 11.1, atol=1e-1)
  @test isapprox(F1, -0.876, atol=1e-3)
  @test isapprox(F2, -0.885, atol=1e-3)
end
