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
lwca = TPT.LWCA(ahs, T, struct=:full)
wca = TPT.WCA(ahs, T)

nfe = TPT.NFE(ρ, zs, TPT.BretonnetSilbert(zs, rc, a))
tb = TPT.WHTB(zd, rd, version=:original)
nfetb = TPT.NFETB(nfe, tb)

sys1 = TPT.TPTSystem(ahs, nfetb, T = 1833., m = 1.018e+5)
sys2 = TPT.TPTSystem(wca, nfetb, m = 1.018e+5)
sys3 = TPT.TPTSystem(lwca, nfetb, m = 1.018e+5)

resid2 = sys2.ref.residue
resid3 = sys2.ref.residue

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
plot([r -> B2(r)*r^2, r -> B3(r)*r^2], 3.0, 5.5, ylims=:auto, labels=["WCA" "LWCA"])

S1 = TPT.structurefactor(sys1.ref)[1,1]
S2 = TPT.structurefactor(sys2.ref)[1,1]
S3 = TPT.structurefactor(sys3.ref)[1,1]
plot([S1, S2, S3], 0.1, 6, labels = ["AHS" "WCA" "LWCA"], ylims=:auto)

# Total free energy
F1 = TPT.helmholtz(sys1)
F2 = TPT.helmholtz(sys2)

@testset "Unit NFETB" begin
  @test isapprox(σ2, 4.24, atol=1e-2)
  @test isapprox(σ3, 4.26, atol=1e-2)
  @test isapprox(F1, -0.876, atol=1e-3)
  @test isapprox(F2, -0.884, atol=1e-3)
end
