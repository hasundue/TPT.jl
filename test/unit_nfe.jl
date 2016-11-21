import TPT

using Base.Test

# using Plots; pyplot()

#
# Single component: Liquid Na
# Ref: I. H. Umar et al.: J. Phys. F: Metal Phys., Vol. 4 (1974), 1691-1706.
#
σ_Na = 6.235

ahs = TPT.AHS(σ = 6.235, ρ = 3.564e-3)
pp = TPT.Ashcroft(1.0, 1.67)
nfe = TPT.NFE(ahs, pp)
sys = TPT.TPTSystem(ahs, nfe, T = 373, m = 4.191e+4)

kF = TPT.fermiwavenumber(nfe)

G = TPT.localfiled(nfe)
# plot(q -> G(q*kF), 0, 6)

ϵ = TPT.dielectric(nfe)
# plot(q -> ϵ(q*kF), 0, 6, ylims=(-1, 10))

ω₀ = TPT.formfactor(nfe)[1]
# plot(q -> ω₀(q*kF), 0, 6, ylims=(-1,1))

F = TPT.wnechar(nfe)[1,1]
# plot(q -> 27.21*F(q*kF), 0, 6, ylims = (-0.1, 0))

u = TPT.pairpotential(nfe)[1,1]
# plot(u, 2, 20, ylims = (-2.1e-3, 2.1e-3))

T = 373.
kB = TPT.kB
U_eg = TPT.internal_eg(nfe)
S_eg = TPT.entropy(nfe, ahs, T)
U_es = TPT.internal_es(nfe, ahs)
U_bs = TPT.internal(ahs, nfe)
K = TPT.kinetic(sys)
S_hs = TPT.entropy(sys) - S_eg
F_hs = -T*kB*S_hs

@testset "Unit NFE" begin
  @test isapprox(kF, 0.473, atol=1e-3)
  @test isapprox(G(kF), 0.267, atol=1e-3)
  @test isapprox(ϵ(kF), 3.46, atol=1e-2)
  @test isapprox(ω₀(kF), -0.141, atol=1e-3)
  @test isapprox(F(kF), -0.0437, atol=1e-4)
  @test isapprox(u(σ_Na), -0.00122, atol=1e-5)

  @test isapprox(U_eg, -0.0816, atol=1e-4)
  @test isapprox(S_eg, 0.0522, atol=1e-4)
  @test isapprox(U_es, -0.147, atol=1e-3)
  @test isapprox(U_bs, -0.00912, atol=1e-5)
  @test isapprox(K, 0.00177, atol=1e-5)
  @test isapprox(F_hs, -0.00856, atol=1e-5)
end
