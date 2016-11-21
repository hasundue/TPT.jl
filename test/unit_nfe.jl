import TPT

using Base.Test

# using Plots; pyplot()

#
# Single component: Liquid Na
# Ref: I. H. Umar et al.: J. Phys. F: Metal Phys., Vol. 4 (1974), 1691-1706.
#
ahs1 = TPT.AHS(σ = 6.087, ρ = 3.564e-3)
pp1 = TPT.Ashcroft(1.0, 1.67)
nfe1 = TPT.NFE(ahs1, pp1)
sys1 = TPT.TPTSystem(ahs1, nfe1, T = 373, m = 4.191e+4)

kF1 = TPT.fermiwavenumber(nfe1)

G1 = TPT.localfiled(nfe1)
# plot(q -> G(q*kF), 0, 6)

ϵ1 = TPT.dielectric(nfe1)
# plot(q -> ϵ(q*kF), 0, 6, ylims=(-1, 10))

ω1₀ = TPT.formfactor(nfe1)[1]
# plot(q -> ω1₀(q*kF), 0, 6, ylims=(-1,1))

F1 = TPT.wnechar(nfe1)[1,1]
# plot(q -> 27.21*F(q*kF), 0, 6, ylims = (-0.1, 0))

u1 = TPT.pairpotential(nfe1)[1,1]
# plot(u1, 2, 20, ylims = (-2.1e-3, 2.1e-3))
u1(6.087)

g = TPT.paircorrelation(sys1)[1,1]
# plot(g, 1, 20)

T = 373.
kB = TPT.kB
U1_eg = TPT.internal_eg(nfe1)
S1_eg = TPT.entropy(nfe1, ahs1, T)
U1_es = TPT.internal_es(nfe1, ahs1)
U1_bs = TPT.internal(ahs1, nfe1)
K1 = TPT.kinetic(sys1)
S1_hs = TPT.entropy(sys1) - S1_eg
F1_hs = -T*kB*S1_hs

#
# Binary
#
σ_Na = 6.235
σ_K = 7.354

ahs2 = TPT.AHS(σ = [σ_Na, σ_K], ρ = 2.461e-3, c = [0.5, 0.5])
pp2 = TPT.Ashcroft([1.0, 1.0], [1.67, 2.12])
nfe2 = TPT.NFE(ahs2, pp2)
sys2 = TPT.TPTSystem(ahs2, nfe2, T = 373, m = [4.191e+4, 4.191e+4])

ω2₀ = TPT.formfactor(nfe2)
# plot(q -> ω2₀(q*kF), 0, 6, ylims=(-1,1))

F2 = TPT.wnechar(nfe2)
# plot(q -> 27.21*F(q*kF), 0, 6, ylims = (-0.1, 0))

u2 = TPT.pairpotential(nfe2)
# plot(u2, 2, 20, ylims = (-2.1e-3, 2.1e-3))


@testset "Unit NFE" begin
  @test isapprox(kF1, 0.473, atol=1e-3)
  @test isapprox(G1(kF), 0.267, atol=1e-3)
  @test isapprox(ϵ1(kF), 3.46, atol=1e-2)
  @test isapprox(ω1₀(kF), -0.141, atol=1e-3)
  @test isapprox(F1(kF), -0.0437, atol=1e-4)
  @test isapprox(u1(6.087), -0.00078, atol=1e-5)

  @test isapprox(U1_eg, -0.0816, atol=1e-4)
  @test isapprox(S1_eg, 0.0522, atol=1e-4)
  @test isapprox(U1_es, -0.147, atol=1e-3)
  @test isapprox(U1_bs, -0.00831, atol=1e-5)
  @test isapprox(K1, 0.00177, atol=1e-5)
  @test isapprox(F1_hs, -0.00921, atol=1e-5)
end
