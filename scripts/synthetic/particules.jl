using Obsea
using NPZ
using Plots
include("../plot.jl")
gr(size = (800, 800))

## Load data & parameters & precompute likelihood
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
tdoalut = TDOALUT(propa, grid)
ℓ = Likelihood(zr, za, tdoalut, models, grid)

## Animation
weights, particles = init(Nt, Np)
@time anim = @animate for t = 100:250
    print(t, "/", Nt, "\r")
    cloud, prevcloud, ℓt = slice(t, particles, ℓ)
    predict!(weights, cloud, prevcloud, ℓt, models, grid)
    update!(weights, cloud, ℓt, models, grid)

    plot()
    scatter!(x[1:t], fill(y, t), grid, color = :blue, markerstrokecolor = :blue)
    plot!(particles, grid, color = :green, alpha = 0.1)
    scatter!([0], [0], label = "OBS", color = :orange)

    @views resample!(weights, particles[1:t, :], ℓt, models, grid)
end
gif(anim, "anim.gif", fps = 10)
