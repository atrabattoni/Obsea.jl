import Obsea:
    parameters,
    TDOALUT,
    likelihood,
    init,
    predict!

using Plots
using Parameters
include("plot.jl")

fs = 50.0
Nfft = 1024
Np = 10000
models, propa, grid = parameters("params.toml", fs, Nfft)
tdoalut = TDOALUT(propa, grid)
zr = rand(grid.Nτ)
zr[160:163] .+= 5.0
zr .^= 2
za = 2π * rand(grid.Nf)
za[101:105] .= 3.0 .+ 2π * rand(5) / 36


ℓ = likelihood(zr, za, tdoalut, models, grid)
plot(grid.r, ℓ.r[:, 1])
heatmap(rad2deg.(grid.a), grid.f, ℓ.a[:, :, 1])

weights, cloud = init(Np)
predict!(weights, cloud, ℓ, models, grid)
plot(cloud, grid)
