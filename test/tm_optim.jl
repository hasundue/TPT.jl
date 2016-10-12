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

Threads.@threads for i in 1:N
    function fopt(rc::Float64)
        nfe = TPT.NFE(p[:ρ][i], p[:T][i], 0.0, p[:zs][i], TPT.Ashcroft(rc))
        nfetb = TPT.NFETB(nfe, tb[i])
        sys = TPT.TPTSystem(wca[i], nfetb)
        g = TPT.prdf(sys)[1,1]

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

@testset "Optimization of rc" begin
  for i in 1:N
    @test Optim.converged(res[i])
    @test 0.9 < rc[i] / rc₀[i] < 1.1
  end
end # testset
