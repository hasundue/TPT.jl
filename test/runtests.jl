using TPT
using Base.Test
using Optim

# some constants
ϵ = eps(Float64)

#
# single-component additive hard-sphere system
#
ahs = AHSSystem(ρ = 1.0, η = 0.45)
σ = ahs.σ[1]
g = prdf(ahs)[1,1]

res = optimize(r -> -g(r), 0.5σ, 1.5σ)
σ_min = Optim.minimizer(res)

@testset "Unary AHS" begin
  @test 0.97σ < σ_min < 1.03σ
end
