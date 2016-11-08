import TPT

using Base.Test

# Auxiliary functions
M2m(M) = M / 1000 / 6.022e+23 / 9.109e-31

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

S = TPT.entropy(ahs)

@testset "AHS" begin
  @test isapprox(S, 9.44, atol=1e-2)
end
