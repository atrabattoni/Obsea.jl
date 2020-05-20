using Obsea
using NPZ
using Plots
include("../plot.jl")
gr()

## Load data & parameters & precompute likelihood
fs = 50.0
Nfft = 1024
Nt = 360
Np = 1000
models, propa, grid = parameters("params.toml", fs, Nfft)
tdoalut = TDOALUT(propa, grid)

zr = randn(grid.Nτ, Nt) .^ 2
za = 2π .* rand(grid.Nf, Nt)

ℓ = Likelihood(zr, za, tdoalut, models, grid)

## Plot raw data
heatmap(1:Nt, grid.r / 1000, log.(ℓ.r[:, :, 1]))
heatmap(
    1:Nt,
    rad2deg.(grid.a),
    log.(sum(ℓ.a, dims = 1)[1, :, :, 1] / size(ℓ.a, 1)),
)
scatter(log.(ℓ.m))
scatter(log.(ℓ.Σm))
scatter(ℓ.m ./ ℓ.Σm)
## Initialize
weights, particles = init(Nt, Np)
t = 0

## Iterate
t += 1
cloud, prevcloud, ℓt = slice(t, particles, ℓ)
plot()
scatter!(particles[t-1, :], grid, label = "prior")

predict!(weights, cloud, prevcloud, ℓt, models, grid)
scatter!(particles[t, :], grid, label = "predict")

update!(weights, cloud, ℓt, models, grid)
scatter!(
    particles[t, :],
    grid,
    label = "update",
    zcolor = weights,
    colorbar = false,
)

@views resample!(weights, particles[1:t, :], ℓ, models, grid)
scatter!(particles[t, :], grid, label = "resample")

q = sum(weights .* (1 .- Obsea.isdead.(cloud)))
plot!(title = string(q))
