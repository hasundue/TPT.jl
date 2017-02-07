import TPT

using Base.Test

using Plots; pyplot()

println("--- Unit WCA ---")

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

B = TPT.blipfunction(sys.ref)[1,1]
plot(B, 4, 8, ylims=:auto)

g_ahs = TPT.paircorrelation(ahs)[1,1]
g_wca = TPT.paircorrelation(sys)[1,1]
σ_wca = TPT.hsdiameter(sys.ref)[1,1]

plot([g_ahs, g_wca], 2, 20, labels = ["AHS" "WCA"], ylims=:auto)

S = TPT.entropy(sys)
U = TPT.internal(sys)
F = TPT.helmholtz(sys)

@testset "Unit WCA" begin
  @test isapprox(σ_wca, 6.17, atol=1e-2)
  @test isapprox(S, 7.56, atol=1e-2)
  @test isapprox(U, -0.236, atol=1e-3)
  @test isapprox(F, -0.244, atol=1e-3)
end
