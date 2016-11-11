import TPT

using Base.Test

# Atomic mass
M = [22.9898, 39.0983]
M2m(M) = M / 1000 / 6.022e+23 / 9.109e-31
m_Na, m_K = map(M2m, M)

#
# Single component: Liquid Na
# Ref: I. H. Umar et al.: J. Phys. F: Metal Phys., Vol. 4 (1974), 1691-1706.
#
ahs1 = TPT.AHSSystem(σ = 6.235, ρ = 3.564e-3)
sys1 = TPT.TPTSystem(ahs1, T = 373, m = m_Na)

K1 = TPT.kinetic(sys1)
S1 = TPT.entropy(sys1)
U1 = TPT.internal(sys1)
F1 = TPT.helmholtz(sys1)

#
# Two components: NaK liquid alloy
# Ref: I. H. Umar et al.: J. Phys. F: Metal Phys., Vol. 4 (1974), 1691-1706.
#
ahs2 = TPT.AHSSystem(ρ = 2.461e-3, σ = [6.235, 7.354], c = [0.5, 0.5])
sys2 = TPT.TPTSystem(ahs2, T = 373, m = [m_Na, m_K])

K2 = TPT.kinetic(sys2)
S2 = TPT.entropy(sys2) # divided by kB
U2 = TPT.internal(sys2)
F2 = TPT.helmholtz(sys2)

@testset "Unit AHS" begin
  @testset "Unary" begin
    @test isapprox(K1, 0.00177, atol=1e-5)
    @test isapprox(S1, 7.25, atol=1e-2)
    @test isapprox(U1, 0.00, atol=1e-2)
    @test isapprox(F1, -0.00679, atol=1e-5)
  end
  @testset "Binary" begin
    @test isapprox(K2, 0.00177, atol=1e-5)
    @test isapprox(S2, 9.44, atol=1e-2)
    @test isapprox(U2, 0.00, atol=1e-2)
    @test isapprox(F2, -0.00938, atol=1e-5)
  end
end
