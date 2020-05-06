using Obsea
using NPZ
using Plots
include("../plot.jl")
gr(size = (800, 800))

## Load data & parameters & precompute likelihood
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
tdoalut = TDOALUT(propa, grid)
ℓr = likelihood(zr, tdoalut, models, grid)
ℓa = likelihood(za, models, grid)
ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)

## Animation
weights, particles = init(Nt, Np)
anim = @animate for t = 1:Nt
    print(t, "/", Nt, "\r")
    cloud, prevcloud, ℓ = slice(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    predict!(weights, cloud, prevcloud, ℓ, models, grid)
    update!(weights, cloud, ℓ, models, grid)

    plot()
    scatter!(x[1:t], fill(y, t), grid, color = :blue, markerstrokecolor = :blue)
    plot!(particles, grid, color = :green, alpha = 0.1)
    scatter!([0], [0], label = "OBS", color = :orange)

    resample!(weights, particles)
end
gif(anim, "anim.gif", fps = 10)
