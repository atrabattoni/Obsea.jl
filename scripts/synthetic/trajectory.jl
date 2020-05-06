using Obsea
using NPZ
using Plots
include("../plot.jl")
plotly()

## Load data & parameters
fs = 50.0
Np = 10000
Nfft = 1024

zr = npzread("ceps_test.npy")
za = npzread("az_test.npy")
x = npzread("x_test.npy")
y = npzread("y_test.npy")
f = npzread("f_test.npy")
Nt = size(zr, 2)

models, propa, grid = parameters("params_test.toml", fs, Nfft)
xe, ye = track(zr, za, models, propa, grid, Np)

## Plot true and estimated track
plot()
scatter!(x, fill(y, Nt), grid; label="true", color=:green)
scatter!(xe, ye, grid; label="estimated", color=:blue)
scatter!([0], [0], label="OBS", color=:orange)
