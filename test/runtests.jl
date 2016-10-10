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

res = optimize(r -> -g(r), 0.9σ, 1.1σ, rel_tol=1e-16)
σ_max = Optim.minimizer(res)

@testset "Unary AHS" begin
  @test σ ≈ σ_max
end
