import TPT
using Base.Test

using DataFrames
import Optim
import NLopt
import Dierckx
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
ahs = Vector{TPT.AHS}(N)
wca = Vector{TPT.LWCA}(N)
# tb = Vector{TPT.WHTB}(N)

for i in 1:N
  ahs[i] = TPT.AHS(σ = p[:σ][i], ρ = p[:ρ][i])
  wca[i] = TPT.LWCA(ahs[i], p[:T][i], struct=:full)
  # tb[i] = TPT.WHTB(p[:zd][i], p[:rd][i], version=:original)
end

a₀ = p[:a]
a = zeros(N)
rd₀ = p[:rd]
rd = zeros(N)
res = zeros(N)
sys = Vector{TPT.TPTSystem}(N)
F = Vector{Vector{Float64}}(N)

#
# Performe optimization of rc
#
# Threads.@threads for i in 1:N
for i in 7
  function foptrd(rd::Float64)::Float64
    @printf "rd = %1.5f\n" rd

    function fopta(a::Float64)::Float64
      local ahs = TPT.AHS(σ = p[:σ][i], ρ = p[:ρ][i])
      local wca = TPT.LWCA(ahs, p[:T][i], struct=:full)
      local pse = TPT.BretonnetSilbert(p[:zs][i], p[:rc][i], a)
      local nfe = TPT.NFE(wca, pse)
      local tb = TPT.WHTB(p[:zd][i], rd, version=:original)
      local nfetb = TPT.NFETB(nfe, tb)
      sys[i] = TPT.TPTSystem(wca, nfetb)
      local g = TPT.paircorrelation(sys[i])[1,1]

      local residue::Float64 = 0
      for j in 1:length(g_exp[i])
        residue += abs(g_exp[i][j][2] - g(g_exp[i][j][1]))
      end
      residue
    end

    opt = Optim.optimize(fopta, 0.25, 0.35, rel_tol = 1e-3)
    a[i] = Optim.minimizer(opt)
    @printf "a = %1.5f\n" a[i]

    F[i] = Vector{Float64}(7)
    Threads.@threads for k in 1:7
      local ρ = (0.7 + (k-1)*0.1) * p[:ρ][i]

      local ahs = TPT.AHS(σ = p[:σ][i], ρ = ρ)
      local wca = TPT.LWCA(ahs, p[:T][i], struct=:full)
      local pse = TPT.BretonnetSilbert(p[:zs][i], p[:rc][i], a[i])
      local nfe = TPT.NFE(wca, pse)
      local tb = TPT.WHTB(p[:zd][i], rd, version=:original)
      local nfetb = TPT.NFETB(nfe, tb)
      sys[i] = TPT.TPTSystem(wca, nfetb, m = p[:m][i])

      F[i][k] = TPT.helmholtz(sys[i])
    end
    @show F[i]

    spl = Dierckx.Spline1D(collect(0.7:0.1:1.3), F[i], k=3, bc="error")
    opt = Optim.optimize(x -> spl(x), 0.7, 1.3)

    residue = abs(1.0 - Optim.minimizer(opt))
    @show residue
  end # foptrd

  fnlopt(x::Vector{Float64}, g)::Float64 = foptrd(x[1])

  # opt = NLopt.Opt(:GN_DIRECT, 1)
  opt = NLopt.Opt(:LN_BOBYQA, 1)
  NLopt.min_objective!(opt, fnlopt)
  NLopt.stopval!(opt, 1e-2)
  NLopt.xtol_abs!(opt, 1e-3)
  NLopt.ftol_abs!(opt, 1e-4)
  NLopt.lower_bounds!(opt, [0.95rd₀[i]])
  NLopt.upper_bounds!(opt, [1.05rd₀[i]])

  (fmin, xmin, result) = NLopt.optimize(opt, [rd₀[i]])

  res[i] = fmin
  rd[i] = xmin[1]
  foptrd(rd[i])
end

#
# Print the numerical results to console
#
println("The optimized values of a:")
for i in 1:N
  @printf "  %2s: %1.5f\n" p[:X][i] a[i]
end

println("The optimized values of rd:")
for i in 1:N
  @printf "  %2s: %1.5f\n" p[:X][i] rd[i]
end

rmin = zeros(N)
println("The positions of minimum of pairpotentials:")
for i in 1:N
  rmin[i] = sys[i].ref.rmin[1,1]
  @printf "  %2s: %1.3f\n" p[:X][i] rmin[i]
end

σ_hs = zeros(N)
println("The effective hard-sphere diameters:")
for i in 7
  σ_hs[i] = sys[i].ref.trial.σ[1]
  @printf "  %2s: %1.3f\n" p[:X][i] σ_hs[i]
end


#
# Save the numerical results as a csv file
#
df = DataFrame(X = p[:X], a = round(a, 5), rd = round(a, 5), res = round(res, 1), rmin = round(rmin, 3), σ_hs = round(σ_hs, 3))
writetable(joinpath(resdir, "results.csv"), df)

# Overwrite the parameters for the subsequent tests
p[:a] = round(a, 5)
p[:rd] = round(a, 5)
writetable(joinpath("data", "parameters", "tm_optim.csv"), p)

#
# Save the graphical results as png files
#
for i in 7
  # 1. Pairpotentials
  u_nfe = TPT.pairpotential(sys[i].pert.nfe)[1,1]
  u_tb = TPT.pairpotential(sys[i].pert.tb)[1,1]
  u_tot = TPT.pairpotential(sys[i].pert)[1,1]

  plot([u_nfe, u_tb, u_tot], 2, 20, ylims=(-0.05, 0.05), labels = ["NFE" "TB" "total"], xlabel="r (a.u.)", ylabel="u(r) (a.u.)")
  vline!([σ_hs[i]], label="HS")
  png(joinpath(resdir, "$(i)-$(p[:X][i])_u"))

  # 2. Pair correlation functions
  g_hs = TPT.paircorrelation(ahs[i])[1,1]
  g_wca = TPT.paircorrelation(sys[i])[1,1]

  plot([g_hs, g_wca], 2, 20, labels=["HS" "WCA"], xlabel="r (a.u.)", ylabel="g(r)", ylims=:auto)
  plot!(g_exp[i], label="exp", xlims=(2,20), ylims=:auto)
  png(joinpath(resdir, "$(i)-$(p[:X][i])_g"))

  # 3. Density dependency of free-energy
  plot(0.7:0.1:1.3, 2625.5*F[i], xlabel="Relative density", ylabel="F (kJ/mol)", ylims=:auto)
  png(joinpath(resdir, "$(i)-$(p[:X][i])_F"))
end


#
# Check the results
#
@testset "TM Optim" begin
  # for i in 1:N
  #   @test isapprox(a[i], a₀[i], atol=1e-3)
  # end
end # testset
