using Obsea
using NPZ
using Plots
gr()

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
tdoalut = TDOALUT(propa, grid)

## Precompute likelihood
ℓr = likelihood(zr, tdoalut, models, grid)
ℓa = likelihood(za, models, grid)
ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)

## Plot raw data
heatmap(collect(1:Nt), grid.τ, zr; color = :viridis)
heatmap(collect(1:Nt), grid.f, za; color = :phase)

## Plot likelihood
heatmap(1:Nt, grid.r / 1000, log.(ℓr[:, :, 1]))
heatmap(1:Nt, rad2deg.(grid.a), log.(sum(ℓa, dims = 1)[1, :, :, 1]))

## Animation
xgrid = -30000.0:500.0:30000.0
cloud = StructArray(State(1, 15.0, x, y, NaN, NaN) for y in xgrid, x in xgrid)
anim = @animate for t = 1:Nt
    print(t, "/", Nt, "\r")
    _, _, ℓ = slice(t, particles, ℓr, ℓa, ℓm, ℓΣm)

    pd = models.pd
    ℓᵣ, ℓₐ = Obsea.make(ℓ, grid)
    r = sqrt.(cloud.x .^ 2 .+ cloud.y .^ 2)
    f = cloud.f
    a = mod.(atan.(cloud.x, cloud.y), 2π)
    m = cloud.m
    w = (1.0 .- pd[m]) .+ pd[m] .* ℓᵣ.(r, m) .* ℓₐ.(f, a, m)

    plot()
    heatmap!(w, xgrid)
    scatter!(x[1:t], fill(y, t), grid, color = :blue, markerstrokecolor = :blue)
    scatter!([0], [0], label = "OBS", color = :orange)
end

gif(anim, "anim.gif", fps = 5)
