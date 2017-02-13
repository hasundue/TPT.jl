import TPT

using Base.Test

using Plots; pyplot()

println("--- Unit LWCA ---")

#
# Lennard-Jones fluid
# Ref: A. Meyer et al. Chemical Physics 49 (1980) 147-152
#
kB = TPT.kB
ϵ = 119.8kB
r₀ = 3.405 / 0.5292
n = 0.844 / r₀^3
T = 0.723 * ϵ / kB

lj = TPT.LennardJones(ϵ, r₀)

u = TPT.pairpotential(lj)[1,1]
rmin = TPT.pairpotential_minimizer(lj)[1,1]
umin = u(rmin)
u′ = TPT.pairpotential_derivative(lj)[1,1]
plot([u, u′], 5, 20, ylims=(umin, -umin))

ahs = TPT.AHS(σ = r₀, ρ = n, approx=:RFA)
lwca = TPT.LWCA(ahs, T, struct=:linear)
wca = TPT.WCA(ahs, T)

sys = TPT.TPTSystem(lwca, lj)
sys2 = TPT.TPTSystem(wca, lj)

σ = sys.ref.trial.σ[1]
σ₀ = 3.474 / 0.5292

σ2 = sys2.ref.trial.σ[1]

B = TPT.blipfunction(sys.ref)[1,1]
B2 = TPT.blipfunction(sys2.ref)[1,1]
plot([B, B2], 5, 8)

g0 = TPT.paircorrelation(ahs)[1,1]
g1 = TPT.paircorrelation(sys)[1,1]
g2 = TPT.paircorrelation(sys2)[1,1]
plot([g0, g1, g2], 5, 15, labels = ["AHS" "LWCA" "WCA"])

S0 = TPT.structurefactor(ahs)[1,1]
S1 = TPT.structurefactor(sys)[1,1]
S2 = TPT.structurefactor(sys2)[1,1]
plot([S0, S1, S2], 0.1, 5, labels=["AHS" "LWCA" "WCA"])

@testset "Unit LWCA" begin
  @test isapprox(σ, 6.55, atol=1e-2)
end
