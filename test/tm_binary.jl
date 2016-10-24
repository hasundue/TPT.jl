import TPT

using Base.Test

using DataFrames
using Optim
using Plots; gr()

# setting up results directory
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_binary")
!isdir(resdir) && mkdir(resdir)

# composition is fixed to 1:1
N = 2
c = [0.5, 0.5]

p = readtable(joinpath("data", "parameters", "parameters.csv"), separator='\t')

data = readtable(joinpath("data", "sf", "Fe-Ni.csv"))

q = data[:q] * 0.5291 # convert Å to a.u.
ndata = length(q)
qmin = q[1]
qmax = q[ndata]

# convert Faber-Ziman to Ashcroft-Langreth:
Sexp = Array{Any}(N,N)
Sexp[1,1] = 1 + √(c[1]*c[1]) * (data[:S11] - 1)
Sexp[1,2] = Sexp[2,1] = √(c[1]*c[2]) * (data[:S12] - 1)
Sexp[2,2] = 1 + √(c[2]*c[2]) * (data[:S22] - 1)

#
# Fitting S with additive hard-sphere
#
function fopt(x::Vector{Float64})::Float64
  ρ = x[1]
  σ = x[2:3]

  ahs = TPT.AHSSystem(ρ = ρ, σ = σ, c = c)
  Sahs = TPT.psf(ahs)

  R = zeros(N,N)

  for i in 1:N, j in 1:N
    i > j && continue
    for k in 1:ndata
      R[i,j] += abs(Sexp[i,j][k] - Sahs[i,j](q[k]))
    end
  end

  for i in 1:N, j in 1:N
    i < j && continue
    R[i,j] = R[j,i]
  end

  return norm(R)
end

σ₀ = [p[:σ][5], p[:σ][7]]
ρ₀ = (p[:ρ][5] + p[:ρ][7]) / 2

x₀ = [ρ₀, σ₀[1], σ₀[2]]
res = optimize(fopt, x₀)

(ρ, σ₁, σ₂) = Optim.minimizer(res)

ahs = TPT.AHSSystem(ρ = ρ, σ = [σ₁, σ₂], c = c)
Sahs = TPT.psf(ahs)

for (i,j) in [(1,1), (1,2), (2,2)]
  plot(q, Sexp[i,j], label="exp")
  plot!(Sahs[i,j], qmin, qmax, label="AHS")
  xlabel!("q (a.u.)")
  ylabel!("S")
  file = string("Fe-Ni_AHS_S", i, j)
  path = joinpath(resdir, file)
  png(path)
end

#
# WCA
#
T = 1769.0

ahs = TPT.AHSSystem(ρ = ρ₀::Float64, σ = σ₀::Vector, c = c::Vector)
pp = TPT.Ashcroft([p[:rc][5], p[:rc][7]])
nfe = TPT.NFE(ρ, c, T, zeros(2), [p[:zs][5], p[:zs][7]], pp)
tb = TPT.WHTB(c, T, 12.0, [p[:zd][5], p[:zd][7]], [p[:rd][5], p[:rd][7]])
nfetb = TPT.NFETB(nfe, tb)

wca = TPT.WCASystem(ahs, T)
sys = TPT.TPTSystem(wca, nfetb)

σ_wca = sys.ref.trial.σ

Swca = TPT.psf(sys)

u_nfe = TPT.pairpotential(sys.pert.nfe)
u_tb = TPT.pairpotential(sys.pert.tb)
u_tot = TPT.pairpotential(sys.pert)

for (i,j) in [(1,1), (1,2), (2,2)]
  # structure factor
  plot(q, Sexp[i,j], label="exp")
  plot!(Swca[i,j], qmin, qmax, label="WCA")
  xlabel!("q (a.u.)")
  ylabel!("S")
  file = string("Fe-Ni_WCA_S", i, j)
  path = joinpath(resdir, file)
  png(path)

  # pair-potential
  plot([u_nfe[i,j], u_tb[i,j], u_tot[i,j]], 2, 20, ylims=(-0.1, 0.1), labels = ["NFE" "TB" "total"], xlabel="r (a.u.)", ylabel="u_$(i)$(j) (r) (a.u.)")
  file = string("Fe-Ni_WCA_u", i, j)
  path = joinpath(resdir, file)
  png(path)
end

@testset "TM Binary" begin
  @testset "AHS" begin
    @test isapprox(ρ, 0.0113, atol=1e-4)
    @test isapprox(σ₁, 4.33, atol=1e-2)
    @test isapprox(σ₂, 4.21, atol=1e-2)
  end
  @testset "AHS-WCA" begin
    @test isapprox(σ_wca[1], 4.53, atol=1e-2)
    @test isapprox(σ_wca[2], 3.80, atol=1e-2)
  end
end
