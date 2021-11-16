# TPT.jl

[![Build Status](https://travis-ci.com/hasundue/TPT.jl.svg?token=Bc5JTjkpgq1XmP3wKNoJ&branch=master)](https://travis-ci.com/hasundue/TPT.jl)

**WARNING: This package was developed just for my personal use and is not maintained anymore.**

TPT.jl is a framework for calculation related to thermodynamic perturbation theory (TPT) in Julia.

This package is at the beginning stage of development, and doesn't have enough features for practical use.

# Example

## Reference system only

Calculate the partial structure factors of a binary additive hard-sphere mixture with the same parameters as Fig. 1 in the classical paper by Ashcroft and Langreth (1967).

```julia
using TPT
system = AHSSystem(η = 0.45, c = [0.8, 0.2], σ = [0.9, 1.0])
S = psf(system)
```

S is a 2x2 matrix of julia functions and you can use your favarite plotting package such as PyPlot and Gadfly to plot them.

```julia
using PyPlot
plot([S[1,1], S[1,2], S[2,2]], 0, 20)
```

## Reference system + Perturbing system (not implemented yet)

```julia
refsys = WCASystem()
pertsys = LJ612System(ϵ = 1.0, σ = 1.0)
system = TPTSystem(refsys, pertsys)
g = prdf(system)
```
