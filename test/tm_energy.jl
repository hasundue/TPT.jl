import TPT

using Base.Test

using DataFrames
using Optim
using Dierckx
using Plots; pyplot()

println("--- TM Energy ---")

# Read elemental parameters (while checking if we are in the correct place)
p = readtable(joinpath("data", "parameters", "tm_optim.csv"))

#
# Prepare the directory to store the results
#
!isdir("results") && mkdir("results")

resdir = joinpath("results", "tm_energy")
!isdir(resdir) && mkdir(resdir)

# Physical constants
kB = TPT.kB

# number of elements available
M = size(p, 1)

#
# expand parameters
#
T = [ float(p[:T][A]) for A in 1:M, B in 1:M ]::Array{Float64,2}

for x in [:σ, :ρ, :m, :zs, :rc, :a, :zd, :rd]
  @eval $x = Array{Vector{Float64},2}(M,M)
  for A in 1:M, B in 1:M
    v = [ p[x][A], p[x][B] ]
    @eval ($x)[$A,$B] = $v
  end
end

# Resulting TPTSystems
res = Array{Vector{TPT.TPTSystem},2}(M,M)

# Various thermodynamic quantities
enr = [:F, :U, :U_nfe, :U_tb, :U_es, :S, :S_gas, :S_conf, :S_ref, :S_nfe, :S_tb]

for E in enr
  ΔE = Symbol(:Δ, E)
  @eval ($E) = Array{Vector{Float64},2}(M,M)
  @eval ($ΔE) = Array{Vector{Float64},2}(M,M)
end

Threads.@threads for B in 1:M
  # (A,B) = ind2sub((M,M), k)

  # a ≥ b && continue

  A = 5

  res[A,B] = Vector{TPT.TPTSystem}(11)

  for E in enr
    ΔE = Symbol(:Δ, E)
    @eval ($E)[$A,$B] = Vector{Float64}(11)
    @eval ($ΔE)[$A,$B] = Vector{Float64}(11)
  end

  for i in 1:11
    x₂ = (i-1) / 10

    ρ₀ = (1-x₂)*ρ[A,B][1] + x₂*ρ[A,B][2]
    c = [1-x₂, x₂]

    ahs = TPT.AHS(ρ = ρ₀, σ = σ[A,B], c = c, approx = "RFA")
    wca = TPT.LWCA(ahs, T[A,B])

    pse = TPT.BretonnetSilbert(zs[A,B], rc[A,B], a[A,B])
    nfe = TPT.NFE(ahs, pse)
    tb = TPT.WHTB(zd[A,B], rd[A,B], c)
    nfetb = TPT.NFETB(nfe, tb)

    sys = TPT.TPTSystem(wca, nfetb, m = m[A,B])

    res[A,B][i] = sys

    U_nfe[A,B][i] = TPT.internal(sys.ref, sys.pert.nfe)
    U_tb[A,B][i] = TPT.internal(sys.ref, sys.pert.tb)
    U_es[A,B][i] = TPT.internal(sys.pert, sys.ref)
    U[A,B][i] = U_nfe[A,B][i] + U_tb[A,B][i] + U_es[A,B][i]

    S_gas[A,B][i] = TPT.entropy_gas(sys)
    S_conf[A,B][i] = TPT.entropy_conf(sys)
    S_ref[A,B][i] = TPT.entropy(sys.ref)
    S_nfe[A,B][i] = TPT.entropy(nfetb.nfe, sys.ref, T[A,B])
    S_tb[A,B][i] = TPT.entropy(nfetb.tb, sys.ref, T[A,B])
    S[A,B][i] = S_gas[A,B][i] + S_conf[A,B][i] + S_ref[A,B][i] + S_nfe[A,B][i]
                + S_tb[A,B][i]

    K = TPT.kinetic(sys) # this can be omitted

    F[A,B][i] = K + U[A,B][i] - kB*T[A,B]*S[A,B][i]
  end

  for i in 1:11
    x₂ = (i-1) / 10
    for E in enr
      ΔE = Symbol(:Δ, E)
      @eval ($ΔE)[$A,$B][$i] =
        ($E)[$A,$B][$i] - (1-$x₂)*($E)[$A,$B][1] - $x₂*($E)[$A,$B][11]
    end
  end
end


#
# Plot various quantities of mixing
#

for A in 5, B in 1:M
  # a ≥ b && continue

  X₁ = p[:X][A]
  X₂ = p[:X][B]
  sys = res[A,B]::Vector{TPT.TPTSystem}

  # settings for plots
  default(xlabel="x($X₂)", ylabel="Energy of mixing (kJ/mol)", ylims=())

  plot() # initialize the plot

  # free energy (total) and internal energy
  for E in [:F, :U, :U_nfe, :U_tb, :U_es]
    ΔE = Symbol(:Δ, E)
    @eval plot!(0:0.1:1.0, 2625.5*($ΔE)[$A,$B], label=$(string(E)))
  end

  # entropic contribution
  for E in [:S, :S_gas, :S_conf, :S_ref, :S_nfe, :S_tb]
    ΔE = Symbol(:Δ, E)
    @eval plot!(0:0.1:1.0, 2625.5*$(T[A,B]*kB)*($ΔE)[$A,$B],
                label=$(string("T", E)))
  end

  # Save the plot as png
  png(joinpath(resdir, "$(X₁)-$(X₂)_F"))

  # Effective hard-sphere diameter
  σ₁ = [ TPT.hsdiameter(sys[i].ref)[1,1] for i in 1:11 ]
  σ₂ = [ TPT.hsdiameter(sys[i].ref)[2,2] for i in 1:11 ]
  σmin = min(minimum(σ₁), minimum(σ₂))
  σmax = max(maximum(σ₁), maximum(σ₂))
  default(ylabel = "Effective HS diameter (a.u.)", ylims=(0.99σmin, 1.01σmax))
  plot(0:0.1:1, σ₁, label=X₁)
  plot!(0:0.1:1, σ₂, label=X₂)
  png(joinpath(resdir, "$(X₁)-$(X₂)_σ"))
end

@testset "TM Energy" begin
end
