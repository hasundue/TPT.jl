import TPT

using Base.Test

# Auxiliary functions
M2m(M) = M / 1000 / 6.022e+23 / 9.109e-31

#
# Liquid Na
# Ref: I. H. Umar et al.: J. Phys. F: Metal Phys., Vol. 4 (1974), 1691-1706.
#
ahs1 = TPT.AHSSystem(1, [6.235], [0.003564], [4.19e+4], 373.)

S1 = TPT.entropy(ahs1)
F1 = TPT.helmholtz(ahs1)

#
# NaK liquid alloy
# Ref: I. H. Umar et al.: J. Phys. F: Metal Phys., Vol. 4 (1974), 1691-1706.
#
N = 2
c = [0.5, 0.5]

σ = [6.235, 7.354]
ρ = 0.002461 * c

M = [22.9898, 39.0983]
m = map(M2m, M)

T = 373.

ahs = TPT.AHSSystem(N, σ, ρ, m, T)

S2 = TPT.entropy(ahs) # divided by kB
F2 = TPT.helmholtz(ahs)

@testset "Unit AHS" begin
  @testset "Unary" begin
    @test isapprox(S1, 7.25, atol=1e-2)
    @test isapprox(F1, -0.00679, atol=1e-5)
  end
  @testset "Binary" begin
    @test isapprox(S2, 9.44, atol=1e-2)
    @test isapprox(F2, -0.00938, atol=1e-5)
  end
end
