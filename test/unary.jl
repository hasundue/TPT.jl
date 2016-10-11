import TPT
using Base.Test

using Optim

ahs = TPT.AHSSystem(ρ = 1.0, η = 0.45)
σ = ahs.σ[1]
g = TPT.prdf(ahs)[1,1]

res = optimize(r -> -g(r), 0.9σ, 1.1σ, rel_tol=1e-16)
σ_max = Optim.minimizer(res)

@testset "Unary AHS" begin
  @test σ ≈ σ_max
end
