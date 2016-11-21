import TPT

using Base.Test

using Plots; pyplot()

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

ahs = TPT.AHS(σ = 4.25, ρ = ρ)
wca = TPT.WCA(ahs, T)

nfe = TPT.NFE(ρ, zs, TPT.BretonnetSilbert(zs, rc, a))
tb = TPT.WHTB(zd, rd)
nfetb = TPT.NFETB(nfe, tb)

sys1 = TPT.TPTSystem(ahs, nfetb, T = 1833., m = 1.018e+5)
sys2 = TPT.TPTSystem(wca, nfetb, m = 1.018e+5)

u_nfe = TPT.pairpotential(nfe)[1,1]
u_tb = TPT.pairpotential(tb)[1,1]
u_tot = TPT.pairpotential(nfetb)[1,1]

plot([u_nfe, u_tb, u_tot], 2, 10, ylims = (-0.06, 0.1), labels = ["NFE" "TB" "Total"])

g1 = TPT.paircorrelation(sys1)[1,1]
g2 = TPT.paircorrelation(sys2)[1,1]
plot([g1, g2], 2, 10, labels = ["AHS" "WCA"])

y2 = TPT.blipfunction(sys2.ref)[1,1]
plot(y2, 3, 6)

U1_ref = TPT.internal(sys1.ref)
U2_ref = TPT.internal(sys2.ref)

U1_nfe = TPT.internal(sys1.ref, sys1.pert.nfe)
U2_nfe = TPT.internal(sys2.ref, sys2.pert.nfe)

U1_tb = TPT.internal(sys1.ref, sys1.pert.tb)
U2_tb = TPT.internal(sys2.ref, sys2.pert.tb)

U1_pair = TPT.internal(sys1.ref, sys1.pert)
U2_pair = TPT.internal(sys2.ref, sys2.pert)

U1_es = TPT.internal(sys1.pert, sys1.ref)
U2_es = TPT.internal(sys2.pert, sys2.ref)

U1 = TPT.internal(sys1)
U2 = TPT.internal(sys2)

S1 = TPT.entropy(sys1)
S2 = TPT.entropy(sys2)

F1 = TPT.helmholtz(sys1)
F2 = TPT.helmholtz(sys2)

@testset "Unit NFETB" begin
  @test isapprox(U1_ref, 0.00, atol=1e-2)
  @test isapprox(U2_ref, 0.0108, atol=1e-4)
  @test isapprox(U1_nfe, 0.0271, atol=1e-4)
  @test isapprox(U2_nfe, 0.0187, atol=1e-4)
  @test isapprox(U1_tb, -0.166, atol=1e-3)
  @test isapprox(U2_tb, -0.171, atol=1e-3)
  @test isapprox(U1_pair, -0.139, atol=1e-3)
  @test isapprox(U2_pair, -0.147, atol=1e-3)
  @test isapprox(U1_es, -0.708, atol=1e-3)
  @test isapprox(U2_es, -0.708, atol=1e-3)
  @test isapprox(U1, -0.847, atol=1e-3)
  @test isapprox(U2, -0.844, atol=1e-3)
  @test isapprox(S1, 10.7, atol=1e-1)
  @test isapprox(S2, 10.8, atol=1e-1)
  @test isapprox(F1, -0.900, atol=1e-3)
  @test isapprox(F2, -0.898, atol=1e-3)
end
