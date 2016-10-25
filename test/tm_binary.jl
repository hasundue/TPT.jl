import TPT

using Base.Test

using DataFrames
using Optim
using Dierckx
using Plots; pyplot()

println("--- TM Binary ---")

#
# setting up the directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_binary")
!isdir(resdir) && mkdir(resdir)

# composition is fixed to 1:1
N = 2
c = [0.5, 0.5]

# Elemental parameters
p = readtable(joinpath("data", "parameters", "tm_optim.csv"), separator='\t')

# number of elements available
M = size(p, 1)

# Arrays to store the results
q_exp = Array{Any,2}(M,M)
S_exp = Array{Any,2}(M,M)
g_exp = Array{Any,2}(M,M)
ahs = Array{TPT.AHSSystem,2}(M,M)
sys = Array{TPT.TPTSystem,2}(M,M)


#
# Processing each binary system
#
Threads.@threads for k in 1:(M^2)
  (a,b) = ind2sub((M,M), k)

  A = p[:X][a]
  B = p[:X][b]

  psffile = joinpath("data", "sf", "$(A)-$(B).csv")

  # Skip if file does not exist
  !isfile(psffile) && continue

  push!(sysname, string(A, '-', B))

  # Initial guess for HS diameters and number density
  T = (p[:T][a] + p[:T][b]) / 2
  σ₀ = [p[:σ][a], p[:σ][b]]
  ρ₀ = (p[:ρ][a] + p[:ρ][b]) / 2

  # Read Sᵢⱼ(q) from csv into dataframe
  data = readtable(psffile)

  # Convert Å to a.u.
  q_exp[a,b] = data[:q] * 0.5291

  ndata = length(q_exp[a,b])
  q_min = q_exp[a,b][1]
  q_max = q_exp[a,b][ndata]

  # Convert Faber-Ziman to Ashcroft-Langreth:
  S_exp[a,b] = Array{Any,2}(N,N)
  S_exp[a,b][1,1] = 1 + √(c[1]*c[1]) * (data[:S11] - 1)
  S_exp[a,b][1,2] = √(c[1]*c[2]) * (data[:S12] - 1)
  S_exp[a,b][2,2] = 1 + √(c[2]*c[2]) * (data[:S22] - 1)

  # We use Faber-Ziman for obtaining RDF
  S = Array{Any,2}(N,N)
  S[1,1] = Spline1D(q_exp[a,b], data[:S11], k=3, bc="zero")
  S[1,2] = Spline1D(q_exp[a,b], data[:S12], k=3, bc="zero")
  S[2,2] = Spline1D(q_exp[a,b], data[:S22], k=3, bc="zero")

  # Convert S(q) to g(r)
  g_exp[a,b] = Array{Function,2}(N,N)
  for (i,j) in [(1,1), (1,2), (2,2)]
    function g(r)::Float64
      val, err = quadgk(q -> (S[i,j](q) - 1) * sin(q*r) / r * q, q_min, q_max)
      1 + val / (2π^2 * ρ₀)
    end
    g_exp[a,b][i,j] = g
  end

  #
  # PART-1. AHS: fitting S with additive hard-sphere
  #
  function fopt(x::Vector{Float64})::Float64
    ρ = x[1]
    σ = x[2:3]

    ahs = TPT.AHSSystem(ρ = ρ, σ = σ, c = c)
    S_ahs = TPT.psf(ahs)

    R::Array{Float64,2} = zeros(N,N)

    for i in 1:N, j in 1:N
      i > j && continue
      for k in 1:ndata
        R[i,j] += abs(S_exp[a,b][i,j][k] - S_ahs[i,j](q_exp[a,b][k]))
      end
    end

    for i in 1:N, j in 1:N
      i < j && continue
      R[i,j] = R[j,i]
    end

    return norm(R)
  end

  x₀ = [ρ₀, σ₀[1], σ₀[2]]
  res = optimize(fopt, x₀)

  (ρ_ahs, σ₁_ahs, σ₂_ahs) = Optim.minimizer(res)
  σ_ahs = [σ₁_ahs, σ₂_ahs]

  ahs[a,b] = TPT.AHSSystem(ρ = ρ_ahs, σ = σ_ahs, c = c::Vector)

  #
  # PART-2. AHS-WCA:
  #
  # Setting up the system
  pp = TPT.Ashcroft([p[:rc][a], p[:rc][b]])
  nfe = TPT.NFE(ρ_ahs, c::Vector, T, zeros(2), [p[:zs][a], p[:zs][b]], pp)
  tb = TPT.WHTB(c::Vector, T, 12.0, [p[:zd][a], p[:zd][b]], [p[:rd][a], p[:rd][b]])
  nfetb = TPT.NFETB(nfe, tb)
  wca = TPT.WCASystem(ahs[a,b], T)

  # Performing WCA optimization
  sys[a,b] = TPT.TPTSystem(wca, nfetb)
end

#
# Prepare DataFrame for tabular output
#
res = DataFrame(System = AbstractString[], ρ_AHS = Float64[], σ₁_AHS = Float64[], σ₂_AHS = Float64[], σ₁_WCA = Float64[], σ₂_WCA = Float64[])


#
# Output the results as graphs while creating DataFrame
#
for a in 1:M, b in 1:M
  A = p[:X][a]
  B = p[:X][b]

  # Skip if file does not exist
  psffile = joinpath("data", "sf", "$(A)-$(B).csv")
  !isfile(psffile) && continue

  println("$(A)-$(B):")

  # Experimental data
  ndata = length(q_exp[a,b])
  q_min = q_exp[a,b][1]
  q_max = q_exp[a,b][ndata]
  #
  # Parameters for AHS
  #
  ρ_ahs = round(sum(ahs[a,b].ρ), 4)
  σ_ahs = round(ahs[a,b].σ, 3)

  println("  AHS:")
  @printf "    %4s: %1.4f\n" "ρ" ρ_ahs
  @printf "    %4s: %1.3f\n" "σ_$A" σ_ahs[1]
  @printf "    %4s: %1.3f\n" "σ_$B" σ_ahs[2]

  #
  # Parameters for AHS-WCA
  #
  σ_wca = round(sys[a,b].ref.trial.σ, 3)

  println("  AHS-WCA:")
  @printf "    %4s: %1.3f\n" "σ_$A" σ_wca[1]
  @printf "    %4s: %1.3f\n" "σ_$B" σ_wca[2]

  # push the results to the DataFrame
  push!(res, ["$(A)-$(B)" ρ_ahs σ_ahs[1] σ_ahs[2] σ_wca[1] σ_wca[2]])

  #
  # Individual uᵢⱼ, Sᵢⱼ and gᵢⱼn comparison with experimental data
  #

  # Sᵢⱼ(q)
  S_ahs = TPT.psf(ahs[a,b])
  S_wca = TPT.psf(sys[a,b])

  # gᵢⱼ(r)
  g_ahs = TPT.prdf(ahs[a,b])
  g_wca = TPT.prdf(sys[a,b])

  # uᵢⱼ(r)
  u_nfe = TPT.pairpotential(sys[a,b].pert.nfe)
  u_tb = TPT.pairpotential(sys[a,b].pert.tb)
  u_tot = TPT.pairpotential(sys[a,b].pert)

  for (i,j) in [(1,1), (1,2), (2,2)]
    # Partial structure factors Sᵢⱼ(q)
    plot(q_exp[a,b], S_exp[a,b][i,j], label="exp")
    plot!(S_ahs[i,j], q_min, q_max, label="AHS")
    plot!(S_wca[i,j], q_min, q_max, label="WCA")
    xlabel!("q (a.u.)")
    ylabel!("S(q)")
    file = string("$(A)-$(B)_S", i, j)
    path = joinpath(resdir, file)
    png(path)

    # Pair distribution functions gᵢⱼ(r)
    plot([g_exp[a,b][i,j], g_ahs[i,j], g_wca[i,j]], 2, 20, labels=["exp" "AHS" "WCA"])
    xlabel!("r (a.u.)")
    ylabel!("g(r)")
    file = string("$(A)-$(B)_g", i, j)
    path = joinpath(resdir, file)
    png(path)

    # Pair-potentials uᵢⱼ(r)
    plot([u_nfe[i,j], u_tb[i,j], u_tot[i,j]], 2, 20, ylims=(-0.1, 0.1), labels = ["NFE" "TB" "total"], xlabel="r (a.u.)", ylabel="u(r) (a.u.)")
    σᵢⱼ = (σ_wca[i] + σ_wca[j]) / 2
    vline!([σᵢⱼ], label="HS")
    file = string("$(A)-$(B)_u", i, j)
    path = joinpath(resdir, file)
    png(path)
  end
end


#
# Output the results as a csv file
#
writetable(joinpath(resdir, "results.csv"), res)


#
# Performing the tests
#

# read the answers
ans = readtable(joinpath("data", "parameters", "tm_binary.csv"))
l = size(ans, 1)

@testset "TM Binary" begin
  @testset "AHS" begin
    for i in 1:l
      @test isapprox(res[:ρ_AHS][i], ans[:ρ_AHS][i], atol=1e-4)
      @test isapprox(res[:σ₁_AHS][i], ans[:σ₁_AHS][i], atol=1e-3)
      @test isapprox(res[:σ₂_AHS][i], ans[:σ₂_AHS][i], atol=1e-3)
    end
  end
  @testset "AHS-WCA" begin
    for i in 1:l
      @test isapprox(res[:σ₁_WCA][i], ans[:σ₁_WCA][i], atol=1e-3)
      @test isapprox(res[:σ₂_WCA][i], ans[:σ₂_WCA][i], atol=1e-3)
    end
  end
end
