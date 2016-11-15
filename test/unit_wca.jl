import TPT

using Base.Test
#
# Single component: Liquid Na
# Ref: I. H. Umar et al.: J. Phys. F: Metal Phys., Vol. 4 (1974), 1691-1706.
#
σ_Na = 6.235

ahs = TPT.AHS(σ = 6.235, ρ = 3.564e-3)
wca = TPT.WCA(ahs, 373)
pseudo = TPT.Ashcroft(1.0, 1.67)
nfe = TPT.NFE(wca, pseudo)
sys = TPT.TPTSystem(wca, nfe, m = 4.191e+4)

g = TPT.paircorrelation(sys)[1,1]
g_σ = g(σ_Na)

# plot(g, 2, 20, labels = ["AHS" "WCA"])

S = TPT.entropy(sys)
U = TPT.internal(sys)
F = TPT.helmholtz(sys)

@testset "Unit WCA" begin
  @test isapprox(g_σ, 2.04, atol=1e-2)
  @test isapprox(S, 7.47, atol=1e-2)
  @test isapprox(U, -0.237, atol=1e-3)
  @test isapprox(F, -0.244, atol=1e-3)
end
