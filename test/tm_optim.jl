import TPT
using Base.Test

using DataFrames
using Optim
using Plots; pyplot()

println("--- TM Optim ---")

# Load elemental parameters
p = readtable(joinpath("data", "parameters", "tm_optim.csv"))
N = size(p, 1) # number of elements

#
# Load experimental g(r) from csv files
#
g_exp = Vector{Vector{Tuple{Float64,Float64}}}(N)
for i in 1:N
    elem = p[:X][i]
    T = p[:T][i]

    data = readtable(joinpath("data", "rdf", "$(elem)$(T).csv"))

    M = size(data, 1)

    g_exp[i] = Vector{Tuple{Float64,Float64}}(M)
    for j in 1:M
        g_exp[i][j] = (data[:r][j] / 0.5291, data[:g][j])
    end
end

p[:T] = convert(DataArrays.DataArray{Float64,1}, p[:T])

#
# Prepare TPTSystem for each system
#
ahs = Vector{TPT.AHS}(N)
wca = Vector{TPT.WCA}(N)
tb = Vector{TPT.WHTB}(N)

for i in 1:N
    ahs[i] = TPT.AHS(σ = p[:σ][i], ρ = p[:ρ][i])
    wca[i] = TPT.WCA(ahs[i], p[:T][i])
    tb[i] = TPT.WHTB(p[:zd][i], p[:rd][i])
end

res = Vector{Any}(N)
rc = Vector{Float64}(N)
sys = Vector{TPT.TPTSystem}(N)


#
# Performe optimization of rc
#
Threads.@threads for i in 1:N
    function fopt(rc::Float64)::Float64
      pse = TPT.Ashcroft(p[:zs][i], rc)
      nfe = TPT.NFE(wca[i], pse)
      nfetb = TPT.NFETB(nfe, tb[i])
      sys[i] = TPT.TPTSystem(wca[i], nfetb)
      g = TPT.paircorrelation(sys[i])[1,1]

      residue::Float64 = 0

      for j in 1:length(g_exp[i])
          residue += abs(g_exp[i][j][2] - g(g_exp[i][j][1]))
      end

      return residue
    end

    res[i] = Optim.optimize(fopt, 1.0, 2.0)
    rc[i] = Optim.minimizer(res[i])
end


#
# Print the numerical results to console
#
println("The optimized values of rc:")
for i in 1:N
  @printf "  %2s: %1.3f\n" p[:X][i] rc[i]
end

rmin = Vector{Float64}(N)
println("The positions of minimum of pairpotentials:")
for i in 1:N
  rmin[i] = sys[i].ref.rmin[1,1]
  @printf "  %2s: %1.3f\n" p[:X][i] rmin[i]
end

σ_hs = Vector{Float64}(N)
println("The effective hard-sphere diameters:")
for i in 1:N
  σ_hs[i] = sys[i].ref.trial.σ[1]
  @printf "  %2s: %1.3f\n" p[:X][i] σ_hs[i]
end


#
# Prepare a directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_optim")
!isdir(resdir) && mkdir(resdir)


#
# Save the numerical results as a csv file
#
df = DataFrame(X = p[:X], rc = round(rc, 3), rmin = round(rmin, 3), σ_hs = round(σ_hs, 3))
writetable(joinpath(resdir, "results.csv"), df)


#
# Save the graphical results as png files
#

# 1. Pairpotentials
for i in 1:N
  u_nfe = TPT.pairpotential(sys[i].pert.nfe)[1,1]
  u_tb = TPT.pairpotential(sys[i].pert.tb)[1,1]
  u_tot = TPT.pairpotential(sys[i].pert)[1,1]

  plot([u_nfe, u_tb, u_tot], 2, 20, ylims=(-0.1, 0.1), labels = ["NFE" "TB" "total"], xlabel="r (a.u.)", ylabel="u(r) (a.u.)")
  vline!([rmin[i]], label="r_min")
  vline!([σ_hs[i]], label="HS dia.")
  png(joinpath(resdir, "$(i)-$(p[:X][i])_u"))
end

# 2. Pair correlation functions
for i in 1:N
  g_hs = TPT.paircorrelation(ahs[i])[1,1]
  g_wca = TPT.paircorrelation(sys[i])[1,1]

  plot([g_hs, g_wca], 2, 20, labels=["HS" "WCA"], xlabel="r (a.u.)", ylabel="g(r)")
  plot!(g_exp[i], label="exp", xlims=(2,20))
  png(joinpath(resdir, "$(i)-$(p[:X][i])_g"))
end


#
# Check the results
#
rc₀ = p[:rc]

@testset "TM Optim" begin
  for i in 1:N
    @test isapprox(rc[i], rc₀[i], atol=2e-3)
  end
end # testset
