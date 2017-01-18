import TPT

using Base.Test

using Plots; pyplot()

println("--- Unit AHS ---")

#
# Prepare a directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "unit_ahs")
!isdir(resdir) && mkdir(resdir)

# default settings for plots
default(xlabel = "r", ylabel = "g(r)", label = "")

#
# A binary hard-sphere mixture
# Ref: S. B. Yuste et al. J. Chem. Phys. 108 (1998) 3683-3693
#
ahs1 = TPT.AHS(η = 0.49, σ = [1.0, 0.3], c = [1/16, 15/16], approx="RFA")

α1 = ahs1.α
g1 = TPT.paircorrelation(ahs1)

# Fig.1(a)
plot(g1[1,1], 1, 1.25, ylims = (0, 10), yticks = 0:1:10)
png(joinpath(resdir, "g11a"))

# Fig.1(b)
plot(g1[1,1], 1, 2.25, xticks = 1:0.2:2.2, ylims = (0.55, 1.25),
     yticks = 0.0:0.1:1.3)
png(joinpath(resdir, "g11b"))

# Fig.2(a)
plot(g1[1,2], 0.65, 0.90, ylims = (0.8, 4.8), yticks = 0:0.5:10)
png(joinpath(resdir, "g12a"))

# Fig.2(b)
plot(g1[1,2], 0.7, 1.9, xticks = (0:0.2:2), ylims = (0.88, 1.18),
     yticks = 0:0.05:2)
png(joinpath(resdir, "g12b"))

# Fig.3(a)
plot(g1[2,2], 0.30, 0.55, xticks = (0.3:0.05:0.6), ylims = (0.8, 3.8), yticks = 0:0.5:4)
png(joinpath(resdir, "g22a"))

# Fig.3(b)
plot(g1[2,2], 0.4, 1.5, xticks = (0:0.2:2), ylims = (0.93, 1.15),
     yticks = 0:0.05:2)
png(joinpath(resdir, "g22b"))

#
# Tests for α values
#
ahs2 = TPT.AHS(η = 0.49, σ = [1.0, 0.3], c = [1/2, 1/2], approx="RFA")
α2 = ahs2.α

ahs3 = TPT.AHS(η = 0.49, σ = [1.0, 0.3, 0.1], c = [1/102, 1/102, 100/102],
               approx="RFA")
α3 = ahs3.α

@testset "Unit AHS" begin
  @test isapprox(α1, 0.01886, atol=1e-5)
  @test isapprox(α2, 0.02784, atol=1e-5)
  @test isapprox(α3, 0.01837, atol=1e-5)
end
