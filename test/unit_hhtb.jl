import TPT

using Base.Test
using Plots; pyplot()

println("--- Unit HHTB ---")

!isdir("results") && mkdir("results")

resdir = joinpath("results", "unit_hhtb")
!isdir(resdir) && mkdir(resdir)

# Parameters (Hartree)
Z = 12
X = ["Ni", "Y"]
V = [10.04, 33.05] / (5.292^3 * 1e-3)
Nd = [8.60, 1.69]
Ed = [-4.93, -1.77] / 27.21
Wd = [3.78, 6.59] / 27.21
r₀ = (V * 3/4π) .^ (1/3)

x₂ = 0.5
x = [1-x₂, x₂]
Ec = mean(Ed)
Ed = Ed - Ec
Ec = 0

Emin = Ec - 0.2
Emax = Ec + 0.2

#
# Single-component case: Ni
#
hhtb1 = TPT.HHTB(Nd[1], Ed[1], Wd[1], r₀[1])

S = TPT.transfermatrix(hhtb1)[1,1]
plot(E -> imag(S(E)), Emin, Emax, xl="E (a.u.)", yl="S")
png(joinpath(resdir, "S_Ni"))

Δ = TPT.selfenergy(hhtb1)[1]
plot(E -> imag(Δ(E)), Emin, Emax, xl="E (a.u.)", yl="Δ")
png(joinpath(resdir, "Δ_Ni"))

Gᵢᵢ = TPT.diagonalgreenfunction(hhtb1)[1,1]
plot(E -> imag(Gᵢᵢ(E)), Emin, Emax, xl="E (a.u.)", yl="Im G_ii")
png(joinpath(resdir, "G_ii_Ni"))

G = TPT.offdiagonalgreenfunction(hhtb1)[1,1]
plot(E -> imag(G(E)), Emin, Emax, xl="E (a.u.)", yl="Im G_ij")
png(joinpath(resdir, "G_ij_Ni"))

D = TPT.densityofstate(hhtb1)
Ef = TPT.fermienergy(hhtb1)
plot(D, Emin, Emax, label="Ni", xl="E (a.u.)", yl="D (states / a.u. atom)")
vline!([Ef], label="Ef")
png(joinpath(resdir, "D_Ni"))

ΔNd = TPT.chargetransfer(hhtb1)[1]

Θ = TPT.bondorder(hhtb1)[1,1]

u = TPT.pairpotential(hhtb1)[1,1]
u′ = TPT.pairpotential_derivative(hhtb1)[1,1]
plot([u, u′], 2, 8, ylims=(-0.1, 0.1), xl="r (a.u.)", yl="u (a.u.)", label=["u" "u′"])
png(joinpath(resdir, "u_Ni"))

E_band = TPT.bandenergy(hhtb1)
E_site = TPT.onsiteenergy(hhtb1)

#
# Binary case: Ni-Y
#
hhtb = TPT.HHTB(x, Nd, Ed, Wd, r₀)

h = TPT.transferintegral(hhtb)
hd = TPT.nntransferintegral(hhtb)

default(xl="E (a.u.)")

S = TPT.transfermatrix(hhtb)
plot(E -> imag(S[1,1](E)), Emin, Emax, label=X[1])
plot!(E -> imag(S[2,2](E)), Emin, Emax, label=X[2])
ylabel!("S")
png(joinpath(resdir, "S"))

Δ = TPT.selfenergy(hhtb)
plot(E -> imag(Δ[1](E)), Emin, Emax)
plot!(E -> imag(Δ[2](E)), Emin, Emax)
ylabel!("Δ")
png(joinpath(resdir, "Δ"))

Gd = TPT.diagonalgreenfunction(hhtb)
plot(E -> imag(Gd[1](E)), Emin, Emax, label=X[1])
plot!(E -> imag(Gd[2](E)), Emin, Emax, label=X[2])
ylabel!("Im G_ii")
png(joinpath(resdir, "G_ii"))

G = TPT.offdiagonalgreenfunction(hhtb)
plot(E -> imag(G[1,1](E)), Emin, Emax, label="$(X[1])-$(X[1])")
plot!(E -> imag(G[1,2](E)), Emin, Emax, label="$(X[1])-$(X[2])")
plot!(E -> imag(G[2,2](E)), Emin, Emax, label="$(X[2])-$(X[2])")
ylabel!("Im G_ij")
png(joinpath(resdir, "G_ij"))

D = TPT.partialdensityofstate(hhtb)
Dtotal = TPT.totaldensityofstate(hhtb)
Ef = TPT.fermienergy(hhtb)
plot(E -> D[1](E), Emin, Emax, label=X[1])
plot!(E -> D[2](E), Emin, Emax, label=X[2])
plot!(E -> Dtotal(E), Emin, Emax, label="total")
vline!([Ef], label="Ef")
ylabel!("D (states / a.u. atom)")
png(joinpath(resdir, "D"))

ΔNd = TPT.chargetransfer(hhtb)

Θ = TPT.bondorder(hhtb)

u = TPT.pairpotential(hhtb)
plot([u[1,1], u[1,2], u[2,2]], 2, 8, labels=["Ni-Ni" "Ni-Y" "Y-Y"], ylims=(-0.1, 0.1), xl="r (a.u.)", yl="u (a.u.)")
png(joinpath(resdir, "u"))

E_band = TPT.bandenergy(hhtb)
E_site = TPT.onsiteenergy(hhtb)

@testset "Unit HHTB" begin
end
