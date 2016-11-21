import TPT

using Base.Test
using Plots; pyplot()

ahs1 = TPT.AHS(η = 0.49, σ = [1.0, 0.3], c = [1/16, 15/16])
g1 = TPT.paircorrelation(ahs1)

default(xlabel = "r", ylabel = "g(r)", label = "")

# Fig.1(a)
g1_1 = TPT.spline(g1[1,1], 1.0, 1.25, 4, bc = "error")
plot(g1_1, 1, 1.25, ylims = (0, 10), yticks = 0:1:10)

# Fig.1(b)
g1_2 = TPT.spline(g1[1,1], 1.0, 2.25, 16, bc = "error")
plot(g1_2, 1, 2.25, xticks = 1:0.2:2.2, ylims = (0.55, 1.25), yticks = 0.0:0.1:1.3)

# Fig.2(a)
# @time plot(g1[1,2], 0.65, 0.90, ylims = (0.8, 4.8), yticks = 0:0.5:10)

# Fig.2(b)
# @time plot(g1[1,2], 0.7, 1.9, xticks = (0:0.2:2), ylims = (0.88, 1.18), yticks = 0:0.05:10)


# ahs3 = TPT.AHS(η = 0.49, σ = [1.0, 0.3, 0.1], c = [1/102, 1/102, 100/102])
# G3 = TPT.paircorrelationlaplace(ahs3)

@testset "Unit AHS" begin
#   @testset "Energy" begin
#     @test isapprox(K1, 0.00177, atol=1e-5)
#     @test isapprox(S1, 7.25, atol=1e-2)
#     @test isapprox(U1, 0.00, atol=1e-2)
#     @test isapprox(F1, -0.00679, atol=1e-5)
#     @test isapprox(K2, 0.00177, atol=1e-5)
#     @test isapprox(S2, 9.44, atol=1e-2)
#     @test isapprox(U2, 0.00, atol=1e-2)
#     @test isapprox(F2, -0.00938, atol=1e-5)
#   end
end
