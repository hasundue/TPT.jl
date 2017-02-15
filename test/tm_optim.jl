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
# Prepare a directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_optim")
!isdir(resdir) && mkdir(resdir)


#
# Load experimental g(r) from csv files
#
g_exp = Vector{Vector{Tuple{Float64,Float64}}}(N)
for i in 1:N
    elem = p[:X][i]
    T = convert(Int, p[:T][i])

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
ahs = [ TPT.AHS(σ = p[:σ][i], ρ = p[:ρ][i]) for i in 1:N ]
wca = [ TPT.LWCA(ahs[i], p[:T][i], struct=:full) for i in 1:N ]
tb = [ TPT.HHTB(p[:zd][i], p[:Ed][i], p[:Wd][i], p[:R0][i]) for i in 1:N ]

rc₀ = p[:rc]
rc = zeros(N)
a₀ = p[:a]
a = zeros(N)
res = zeros(N)
sys = Vector{TPT.TPTSystem}(N)
F = Vector{Vector{Float64}}(N)

#
# Performe optimization of rc
#
Threads.@threads for i in 1:N
  function fopt(rc::Float64)::Float64
    pse = TPT.Ashcroft(p[:zs][i], rc)
    nfe = TPT.NFE(wca[i], pse, :IU)
    nfetb = TPT.NFETB(nfe, tb[i])
    sys[i] = TPT.TPTSystem(wca[i], nfetb)

    g = TPT.paircorrelation(sys[i])[1,1]

    residue::Float64 = 0

    for j in 1:length(g_exp[i])
      residue += abs(g_exp[i][j][2] - g(g_exp[i][j][1]))
    end

    return residue
  end

  # fnlopt(x::Vector{Float64}, g)::Float64 = fopt(x[1])
  #
  # opt = NLopt.Opt(:GN_DIRECT, 1)
  # NLopt.min_objective!(opt, fnlopt)
  # NLopt.stopval!(opt, 1e-0)
  # NLopt.xtol_abs!(opt, 1e-3)
  # NLopt.ftol_abs!(opt, 1e-4)
  # NLopt.lower_bounds!(opt, [0.5rc₀[i]])
  # NLopt.upper_bounds!(opt, [1.5rc₀[i]])
  #
  # (fmin, xmin, result) = NLopt.optimize(opt, [rc₀[i]])
  #
  # res[i] = result
  # rc[i] = xmin[1]
  # fopt(xmin)

  result = Optim.optimize(fopt, 0.5, 2.0, rel_tol = 1e-3)
  rc[i] = Optim.minimizer(result)
  res[i] = Optim.minimum(result)

  fopt(rc[i])

  #
  # Density dependency of free-energy
  #
  # F[i] = Vector{Float64}(5)
  # for k in 1:5
  #   ρ = (0.8 + (k-1)*0.1) * p[:ρ][i]
  #
  #   ahs = TPT.AHS(σ = p[:σ][i], ρ = ρ)
  #   wca = TPT.LWCA(ahs, p[:T][i], struct=:full)
  #   tb = TPT.WHTB(p[:zd][i], p[:rd][i], version=:original)
  #   pse = TPT.BretonnetSilbert(p[:zs][i], p[:rc][i], a[i])
  #   nfe = TPT.NFE(wca, pse)
  #   nfetb = TPT.NFETB(nfe, tb)
  #   sys = TPT.TPTSystem(wca, nfetb, m = p[:m][i])
  #
  #   F[i][k] = TPT.helmholtz(sys)
  # end
end


#
# Print the numerical results to console
#
println("The optimized values of rc:")
for i in 1:N
  @printf "  %2s: %1.3f\n" p[:X][i] rc[i]
end

rmin = zeros(N)
println("The positions of minimum of pairpotentials:")
for i in 1:N
  rmin[i] = sys[i].ref.rmin[1,1]
  @printf "  %2s: %1.3f\n" p[:X][i] rmin[i]
end

σ_hs = zeros(N)
println("The effective hard-sphere diameters:")
for i in 1:N
  σ_hs[i] = sys[i].ref.trial.σ[1]
  @printf "  %2s: %1.3f\n" p[:X][i] σ_hs[i]
end


#
# Save the numerical results as a csv file
#
df = DataFrame(X = p[:X], rc = round(a, 3), res = round(res, 1), rmin = round(rmin, 3), σ_hs = round(σ_hs, 3))
writetable(joinpath(resdir, "results.csv"), df)

# Overwrite the parameters for the subsequent tests
p[:rc] = round(rc, 3)
writetable(joinpath("data", "parameters", "tm_optim.csv"), p)

#
# Save the graphical results as png files
#

# 1. Pairpotentials
for i in 1:N
  u_nfe = TPT.pairpotential(sys[i].pert.nfe)[1,1]
  u_tb = TPT.pairpotential(sys[i].pert.tb)[1,1]
  u_tot = TPT.pairpotential(sys[i].pert)[1,1]

  plot([u_nfe, u_tb, u_tot], 2, 20, ylims=(-0.05, 0.05), labels = ["NFE" "TB" "total"], xlabel="r (a.u.)", ylabel="u(r) (a.u.)")
  vline!([σ_hs[i]], label="HS")
  png(joinpath(resdir, "$(i)-$(p[:X][i])_u"))
end

# 2. Pair correlation functions
for i in 1:N
  g_hs = TPT.paircorrelation(ahs[i])[1,1]
  g_wca = TPT.paircorrelation(sys[i])[1,1]

  plot([g_hs, g_wca], 2, 20, labels=["HS" "WCA"], xlabel="r (a.u.)", ylabel="g(r)", ylims=:auto)
  plot!(g_exp[i], label="exp", xlims=(2,20), ylims=:auto)
  png(joinpath(resdir, "$(i)-$(p[:X][i])_g"))
end

# 3. Density dependency of free-energy
# for i in 1:N
#   plot(0.8:0.1:1.2, 2625.5*F[i], xlabel="Relative density", ylabel="F (kJ/mol)", ylims=:auto)
#   png(joinpath(resdir, "$(i)-$(p[:X][i])_F"))
# end


#
# Check the results
#
@testset "TM Optim" begin
  # for i in 1:N
  #   @test isapprox(a[i], a₀[i], atol=1e-3)
  # end
end # testset
