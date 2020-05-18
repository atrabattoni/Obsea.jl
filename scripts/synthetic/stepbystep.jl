using Obsea
using NPZ
using Plots
include("../plot.jl")
plotly()

## Load data & parameters & precompute likelihood
fs = 50.0
Np = 1000
Nfft = 1024

zr = npzread("data/ceps.npy")
za = npzread("data/az.npy")
x = npzread("data/x.npy")
y = npzread("data/y.npy")
f = npzread("data/f.npy")
Nt = size(zr, 2)

models, propa, grid = parameters("params.toml", fs, Nfft)
tdoalut = TDOALUT(propa, grid)

ℓr = likelihood(zr, tdoalut, models, grid)
ℓa = likelihood(za, models, grid)
ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)
ℓ = Obsea.Likelihood(ℓr, ℓa, ℓm, ℓΣm)

## Initialize
weights, particles = init(Nt, Np)
t = 30
## Iterate
t += 1
cloud, prevcloud, ℓt = slice(t, particles, ℓ)
plot()
scatter!(particles[t-1, :], grid, label="prior")

predict!(weights, cloud, prevcloud, ℓt, models, grid)
scatter!(particles[t, :], grid, label="predict")

update!(weights, cloud, ℓt, models, grid)
scatter!(particles[t, :], grid, label="update", zcolor=weights, colorbar=false)

@views resample!(weights, particles[1:t, :], ℓ, models, grid)
scatter!(particles[t, :], grid, label="resample")
