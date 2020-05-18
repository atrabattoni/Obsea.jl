using Obsea
using NPZ
using Plots
include("../plot.jl")
plotly()

## Load data & parameters
fs = 50.0
Np = 10000
Nfft = 1024

zr = npzread("data/ceps.npy")
za = npzread("data/az.npy")
x = npzread("data/x.npy")
y = npzread("data/y.npy")
f = npzread("data/f.npy")
Nt = size(zr, 2)

models, propa, grid = parameters("params.toml", fs, Nfft)
xe, ye = track(zr, za, models, propa, grid, Np)

## Plot true and estimated track
plot()
scatter!(x, fill(y, Nt), grid; label="true", color=:green)
scatter!(xe, ye, grid; label="estimated", color=:blue)
scatter!([0], [0], label="OBS", color=:orange)
