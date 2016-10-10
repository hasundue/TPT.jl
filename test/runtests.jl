using TPT
using Base.Test

# some constants
ϵ = eps(Float64)

#
# single-component additive hard-sphere system
#
ahs = AHSSystem(ρ = 1.0, η = 0.45)
σ = ahs.σ[1]
g = prdf(ahs)[1,1]

@testset "Unary AHS" begin
  @test g(σ) ≈ g(σ + ϵ)
  @test g(σ - ϵ) ≈ 0.0
end
