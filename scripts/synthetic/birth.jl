import Obsea:
    parameters,
    TDOALUT,
    likelihood,
    marginalize,
    Likelihood,
    init,
    predict!,
    slice

using Plots
using Parameters
include("../plot.jl")

fs = 50.0
Nfft = 1024
Np = 10000
models, propa, grid = parameters("params.toml", fs, Nfft)
tdoalut = TDOALUT(propa, grid)
zr = rand(grid.Nτ, 1)
zr[160:163] .+= 5.0
zr .^= 2
za = 2π * rand(grid.Nf, 1)
za[101:105] .= 3.0 .+ 2π * rand(5) / 36


tdoalut = TDOALUT(propa, grid)
ℓ = Likelihood(zr, za, tdoalut, models, grid)

plot(grid.r, ℓ.r[:, 1])
heatmap(rad2deg.(grid.a), grid.f, ℓ.a[:, :, 1])

weights, particles = init(Np, 1)
cloud, prevcloud, ℓt = slice(1, particles, ℓ)
predict!(weights, cloud, prevcloud, ℓt, models, grid)
plot()
scatter!(cloud, grid)
