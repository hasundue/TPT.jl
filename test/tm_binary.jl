import TPT

using Base.Test

using DataFrames
using Optim
using Dierckx
using Plots; pyplot()

println("--- TM Binary ---")

# Elemental parameters
p = readtable(joinpath("data", "parameters", "tm_optim.csv"))

# Binary parameters
ans = readtable(joinpath("data", "parameters", "tm_binary.csv"))

# number of elements available
M = size(p, 1)

#
# setting up the directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_binary")
!isdir(resdir) && mkdir(resdir)

# composition is fixed to 1:1
N = 2
c = [0.5, 0.5]

# Arrays to store the results
isdata = falses(M,M)
data = Array{Any,2}(M,M)
q_exp = Array{Any,2}(M,M)
S_exp = Array{Any,2}(M,M)
g_exp = Array{Any,2}(M,M)
ahs = Array{TPT.AHS,2}(M,M)
sys = Array{TPT.TPTSystem,2}(M,M)


#
# Load experimental data
#
for a in 1:M, b in 1:M
  A = p[:X][a]
  B = p[:X][b]

  psffile = joinpath("data", "sf", "$(A)-$(B).csv")

  # Skip if file does not exist
  !isfile(psffile) && continue

  # Remember that this system has experimental data
  isdata[a,b] = true

  # Read Sᵢⱼ(q) from csv into dataframe
  data[a,b] = readtable(psffile)

  # Convert Å to a.u.
  q_exp[a,b] = data[a,b][:q] * 0.5291
end

#
# Processing each binary system
#
Threads.@threads for k in 1:(M^2)
  (a,b) = ind2sub((M,M), k)

  # skip if experimental data does not exist
  !isdata[a,b] && continue

  # names of components
  A = p[:X][a]
  B = p[:X][b]

  # infomation about experimental data
  ndata = length(q_exp[a,b])
  q_min = q_exp[a,b][1]
  q_max = q_exp[a,b][ndata]

  # Temperature
  entry = ans[ans[:System] .== "$(A)-$(B)", :T]
  T = entry[1]

  # Initial guess for HS diameters and number density
  σ₀ = [p[:σ][a], p[:σ][b]]
  ρ₀ = (p[:ρ][a] + p[:ρ][b]) / 2

  # Convert Faber-Ziman to Ashcroft-Langreth:
  S_exp[a,b] = Array{Any,2}(N,N)
  S_exp[a,b][1,1] = 1 + √(c[1]*c[1]) * (data[a,b][:S11] - 1)
  S_exp[a,b][1,2] = √(c[1]*c[2]) * (data[a,b][:S12] - 1)
  S_exp[a,b][2,2] = 1 + √(c[2]*c[2]) * (data[a,b][:S22] - 1)

  # We use Faber-Ziman for obtaining RDF
  S = Array{Any,2}(N,N)
  S[1,1] = Spline1D(q_exp[a,b], data[a,b][:S11], k=3, bc="zero")
  S[1,2] = Spline1D(q_exp[a,b], data[a,b][:S12], k=3, bc="zero")
  S[2,2] = Spline1D(q_exp[a,b], data[a,b][:S22], k=3, bc="zero")

  # Convert S(q) to g(r)
  g_exp[a,b] = Array{Function,2}(N,N)
  for (i,j) in [(1,1), (1,2), (2,2)]
    g_exp[a,b][i,j] = r -> begin
      val, err = quadgk(q -> (S[i,j](q) - 1) * sin(q*r) / r * q, q_min, q_max)
      1 + val / (2π^2 * ρ₀)
    end
  end

  #
  # PART-1. AHS: fitting S with additive hard-sphere
  #
  function fopt(σ::Vector{Float64})::Float64
    ahs = TPT.AHS(ρ = ρ₀, σ = σ, c = c)
    S_ahs = TPT.structurefactor(ahs)

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

  opt = Optim.optimize(fopt, σ₀, g_tol = 1e-3)

  (σ₁_ahs, σ₂_ahs) = Optim.minimizer(opt)
  σ_ahs = [σ₁_ahs, σ₂_ahs]

  ahs[a,b] = TPT.AHS(ρ = ρ₀, σ = σ_ahs, c = c::Vector)

  #
  # PART-2. AHS-WCA:
  #
  m = [p[:m][a], p[:m][b]]
  zs = [p[:zs][a], p[:zs][b]]
  rc = [p[:rc][a], p[:rc][b]]
  pa = [p[:a][a], p[:a][b]]
  zd = [p[:zd][a], p[:zd][b]]
  rd = [p[:rd][a], p[:rd][b]]

  # pse = TPT.Ashcroft(zs, rc)
  pse = TPT.BretonnetSilbert(zs, rc, pa)
  nfe = TPT.NFE(ahs[a,b], pse)
  tb = TPT.WHTB(zd, rd, c, version=:modified)
  nfetb = TPT.NFETB(nfe, tb)
  # wca = TPT.LWCA(ahs[a,b], T, struct=:full)
  wca = TPT.WCA(ahs[a,b], T)

  # Performing WCA optimization
  sys[a,b] = TPT.TPTSystem(wca, nfetb, m = m)
end

#
# Prepare DataFrame for tabular output
#
res = DataFrame(System = AbstractString[], T = Float64[], ρ_ahs = Float64[], σ₁_ahs = Float64[], σ₂_ahs = Float64[], σ₁_wca = Float64[], σ₂_wca = Float64[])


#
# Output the results as graphs while creating DataFrame
#
for a in 1:M, b in 1:M
  A = p[:X][a]
  B = p[:X][b]

  # Skip if data does not exist
  !isdata[a,b] && continue

  println("$(A)-$(B):")

  # Temperature
  T = ans[ans[:System] .== "$(A)-$(B)", :T][1]
  @printf "  T: %4.0f K\n" T

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
  @printf "    ρ: %1.4f\n" ρ_ahs
  @printf "    %-4s: %1.3f\n" "σ_$A" σ_ahs[1]
  @printf "    %-4s: %1.3f\n" "σ_$B" σ_ahs[2]

  #
  # Parameters for AHS-WCA
  #
  σ_wca = round(sys[a,b].ref.trial.σ, 3)

  println("  AHS-WCA:")
  @printf "    %-4s: %1.3f\n" "σ_$A" σ_wca[1]
  @printf "    %-4s: %1.3f\n" "σ_$B" σ_wca[2]

  # push the results to the DataFrame
  push!(res, ["$(A)-$(B)" T ρ_ahs σ_ahs[1] σ_ahs[2] σ_wca[1] σ_wca[2]])

  #
  # Individual uᵢⱼ, Sᵢⱼ and gᵢⱼn comparison with experimental data
  #

  # Sᵢⱼ(q)
  S_ahs = TPT.structurefactor(ahs[a,b])
  S_wca = TPT.structurefactor(sys[a,b])

  # gᵢⱼ(r)
  g_ahs = TPT.paircorrelation(ahs[a,b])
  g_wca = TPT.paircorrelation(sys[a,b])

  # uᵢⱼ(r)
  u_nfe = TPT.pairpotential(sys[a,b].pert.nfe)
  u_tb = TPT.pairpotential(sys[a,b].pert.tb)
  u_tot = TPT.pairpotential(sys[a,b].pert)

  # Bᵢⱼ(r)
  Bl = TPT.blipfunction(sys[a,b].ref)

  for (i,j) in [(1,1), (1,2), (2,2)]
    # Partial structure factors Sᵢⱼ(q)
    plot(q_exp[a,b], S_exp[a,b][i,j], label="exp", ylims=:auto)
    plot!(S_ahs[i,j], q_min, q_max, label="AHS")
    plot!(S_wca[i,j], q_min, q_max, label="WCA")
    xlabel!("q (a.u.)")
    ylabel!("S(q)")
    file = string("$(a)-$(b)_$(A)-$(B)_S", i, j)
    path = joinpath(resdir, file)
    png(path)

    # Pair distribution functions gᵢⱼ(r)
    plot([g_exp[a,b][i,j], g_ahs[i,j], g_wca[i,j]], 2, 20, labels=["exp" "AHS" "WCA"], ylims=(0, 4.5))
    xlabel!("r (a.u.)")
    ylabel!("g(r)")
    file = string("$(a)-$(b)_$(A)-$(B)_g", i, j)
    path = joinpath(resdir, file)
    png(path)

    # Pair-potentials uᵢⱼ(r)
    plot([u_nfe[i,j], u_tb[i,j], u_tot[i,j]], 2, 20, ylims=(-0.1, 0.1), labels = ["NFE" "TB" "total"], xlabel="r (a.u.)", ylabel="u(r) (a.u.)")
    σᵢⱼ = (σ_wca[i] + σ_wca[j]) / 2
    vline!([σᵢⱼ], label="HS")
    file = string("$(a)-$(b)_$(A)-$(B)_u", i, j)
    path = joinpath(resdir, file)
    png(path)

    # blip function B(r)
    plot(Bl[i,j], 3, 6, xlabel="r (a.u.)", ylabel="B(r)", label = "", ylims=:auto)
    file = string("$(a)-$(b)_$(A)-$(B)_B", i, j)
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
l = size(ans, 1)

@testset "TM Binary" begin
  # @testset "AHS" begin
  #   for i in 1:l
  #     @test isapprox(res[:ρ_ahs][i], ans[:ρ_ahs][i], atol=1e-4)
  #     @test isapprox(res[:σ₁_ahs][i], ans[:σ₁_ahs][i], atol=1e-3)
  #     @test isapprox(res[:σ₂_ahs][i], ans[:σ₂_ahs][i], atol=1e-3)
  #   end
  # end
  # @testset "AHS-WCA" begin
  #   for i in 1:l
  #     @test isapprox(res[:σ₁_wca][i], ans[:σ₁_wca][i], atol=1e-3)
  #     @test isapprox(res[:σ₂_wca][i], ans[:σ₂_wca][i], atol=1e-3)
  #   end
  # end
end
