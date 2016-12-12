import TPT

using Base.Test

using DataFrames
using Optim
using Dierckx
using Plots; pyplot()

println("--- TM Energy ---")

# Read elemental parameters
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

# Resulting TPTSystems
res = Array{Vector{TPT.TPTSystem},2}(M,M)

# Various thermodynamic quantities
enr = [:F, :U, :U_nfe, :U_tb, :U_es, :S, :S_gas, :S_conf, :S_ref, :S_nfe, :S_tb]

for E in enr
  ΔE = Symbol(:Δ, E)
  @eval $(ΔE) = Array{Vector{Float64},2}(M,M)
end

Threads.@threads for b in 1:M
  # (A,B) = ind2sub((M,M), k)

  # a ≥ b && continue

  A = 5

  # Temperature
  T::Float64 = p[:T][A]

  # The other parameters
  for x in [:σ, :ρ, :m, :zs, :rc, :a, :zd, :rd]
    @eval $x = [ $(p[x][A]), $(p[x][B]) ]
  end

  res[A,B] = Vector{TPT.TPTSystem}(11)

  for E in enr
    ΔE = Symbol(:Δ, E)
    @eval $(ΔE)[$A,$B] = Vector{Float64}(11)
    @eval $E = Vector{Float64}(11)
  end

  for i in 1:11
    x₂ = (i-1) / 10

    ρ₀ = (1-x₂)*ρ[1] + x₂*ρ[2]

    ahs = TPT.AHS(ρ = ρ₀, σ = σ, c = [1-x₂, x₂], approx = "RFA")
    wca = TPT.LWCA(ahs, T)

    pse = TPT.BretonnetSilbert(zs, rc, a)
    nfe = TPT.NFE(ahs, pse)
    tb = TPT.WHTB(zd, rd)
    nfetb = TPT.NFETB(nfe, tb)

    sys = TPT.TPTSystem(wca, nfetb, m = m)

    res[A,B][i] = sys

    U_nfe[i] = TPT.internal(sys.ref, sys.pert.nfe)
    U_tb[i] = TPT.internal(sys.ref, sys.pert.tb)
    U_es[i] = TPT.internal(sys.pert, sys.ref)
    U[i] = U_nfe[i] + U_tb[i] + U_es[i]

    S_gas[i] = TPT.entropy_gas(sys)
    S_conf[i] = TPT.entropy_conf(sys)
    S_ref[i] = TPT.entropy(sys.ref)
    S_nfe[i] = TPT.entropy(nfetb.nfe, sys.ref, T)
    S_tb[i] = TPT.entropy(nfetb.tb, sys.ref, T)
    S[i] = S_gas[i] + S_conf[i] + S_ref[i] + S_nfe[i] + S_tb[i]

    K = TPT.kinetic(sys) # this can be omitted

    F[i] = K + U[i] - kB*T*S[i]
  end

  for i in 1:11
    x₂ = (i-1) / 10
    for E in enr
      ΔE = Symbol(:Δ, E)
      @eval ($ΔE)[$A,$B][$i] = $(E)[$i] - (1-$x₂)*($E)[1] - $x₂*($E)[11]
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

  T::Float64 = p[:T][A]

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
    @eval plot!(0:0.1:1.0, 2625.5*$(T*kB)*($ΔE)[$A,$B],
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
