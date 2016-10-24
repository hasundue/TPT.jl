import TPT
using Base.Test

using DataFrames
using Optim

p = readtable(joinpath("data", "parameters", "parameters.csv"), separator='\t')

N = size(p, 1) # number of elements

g_exp = Vector{Vector{Tuple{Float64,Float64}}}(N)

for i in 1:N
    elem = p[:X][i]
    T = p[:T][i]

    data = readtable(joinpath("data", "rdf", "$(elem)$(T).csv"), separator = '\t')

    M = size(data, 1)

    g_exp[i] = Vector{Tuple{Float64,Float64}}(M)
    for j in 1:M
        g_exp[i][j] = (data[:r][j] / 0.5291, data[:g][j])
    end
end

p[:T] = convert(DataArrays.DataArray{Float64,1}, p[:T])

ahs = Vector{TPT.AHSSystem}(N)
wca = Vector{TPT.WCASystem}(N)
tb = Vector{TPT.WHTB}(N)

for i in 1:N
    ahs[i] = TPT.AHSSystem(σ = p[:σ][i], ρ = p[:ρ][i])
    wca[i] = TPT.WCASystem(ahs[i], p[:T][i])
    tb[i] = TPT.WHTB(p[:T][i], p[:zd][i], p[:rd][i])
end

res = Vector{Any}(N)
rc = Vector{Float64}(N)
sys = Vector{TPT.TPTSystem}(N)

Threads.@threads for i in 1:N
    function fopt(rc::Float64)
        nfe = TPT.NFE(p[:ρ][i], p[:T][i], 0.0, p[:zs][i], TPT.Ashcroft(rc))
        nfetb = TPT.NFETB(nfe, tb[i])
        sys[i] = TPT.TPTSystem(wca[i], nfetb)
        g = TPT.prdf(sys[i])[1,1]

        residue::Float64 = 0.0

        for j in 1:length(g_exp[i])
            residue += abs(g_exp[i][j][2] - g(g_exp[i][j][1]))
        end

        return residue
    end

    res[i] = optimize(fopt, 1.0, 2.0)
    rc[i] = Optim.minimizer(res[i])
end

rc₀ = p[:rc]

println("The optimized values of rc:")
for i in 1:N
  @printf "  %2s: %1.3f\n" p[:X][i] rc[i]
end

println("The positions of minimum of pairpotentials:")
for i in 1:N
  @printf "  %2s: %1.3f\n" p[:X][i] sys[i].ref.rmin[1,1]
end

println("The effective hard-sphere diameters:")
for i in 1:N
  @printf "  %2s: %1.3f\n" p[:X][i] sys[i].ref.trial.σ[1]
end

resdir = joinpath("results", "tm_optim", "u")
!isdir(resdir) && mkdir(resdir)

for i in 1:N
  # pairpotentials
  u_nfe = TPT.pairpotential(sys[i].pert.nfe)[1,1]
  u_tb = TPT.pairpotential(sys[i].pert.tb)[1,1]
  u_tot = TPT.pairpotential(sys[i].pert)[1,1]

  rmin::Float64 = sys[i].ref.rmin[1,1]
  σhs::Float64 = sys[i].ref.trial.σ[1]

  plot([u_nfe, u_tb, u_tot], 2, 20, ylims=(-0.1, 0.1), labels = ["NFE" "TB" "total"], xlabel="r (a.u.)", ylabel="u(r) (a.u.)")
  vline!([rmin], label="r_min")
  vline!([σhs], label="HS dia.")
  png(joinpath(resdir, "$(i)-$(p[:X][i])"))
end

resdir = joinpath("results", "tm_optim", "g")
!isdir(resdir) && mkdir(resdir)

for i in 1:N
  g_hs = TPT.prdf(ahs[i])[1,1]
  g_wca = TPT.prdf(sys[i])[1,1]

  plot([g_hs, g_wca], 2, 20, labels=["HS" "WCA"], xlabel="r (a.u.)", ylabel="g(r)")
  plot!(g_exp[i], label="exp", xlims=(2,20))
  png(joinpath(resdir, "$(i)-$(p[:X][i])"))
end

@testset "TM Optim" begin
  for i in 1:N
    @test isapprox(rc[i], rc₀[i], atol=1e-3)
  end
end # testset
