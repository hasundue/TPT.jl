import TPT

using Base.Test

#
# Lennard-Jones
# ref: A. Meyer et al. CHemical Physics 49 (1980) 147-152
#
kB = TPT.kB
ϵ = 119.8kB
r₀ = 3.405 / 0.5292
n = 0.844 / r₀^3
T = 0.723 * ϵ / kB

lj = TPT.LennardJones(ϵ, r₀)

u = TPT.pairpotential(lj)[1,1]
plot(u, 5, 20, ylims=(-1e-3, 1e-3))

ahs = TPT.AHS(σ = r₀, ρ = n, approx="RFA")
lwca = TPT.LWCA(ahs, 86.6)

sys = TPT.TPTSystem(lwca, lj)

σ = sys.ref.trial.σ[1]
σ₀ = 3.474 / 0.5292

B = TPT.blipfunction(sys.ref)[1,1]
plot(r -> B(r)*r^2, 6, 7)

g0 = TPT.paircorrelation(ahs)[1,1]
g1 = TPT.paircorrelation(sys)[1,1]

plot([g0, g1], 5, 20)

@testset "Unit LWCA" begin
  @test isapprox(σ, 6.55, atol=1e-2)
end
