using Obsea
using NPZ
using Plots
using StructArrays
include("../plot.jl")
gr()

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
tdoalut = TDOALUT(propa, grid)

## Precompute likelihood
ℓ = Likelihood(zr, za, tdoalut, models, grid)

## Plot raw data
heatmap(collect(1:Nt), grid.τ, zr; color = :viridis)
heatmap(collect(1:Nt), grid.f, za; color = :phase)

## Plot likelihood
heatmap(1:Nt, grid.r / 1000, log.(ℓr[:, :, 1]))
heatmap(1:Nt, rad2deg.(grid.a), log.(sum(ℓa, dims = 1)[1, :, :, 1]))

## Animation
xgrid = -30000.0:500.0:30000.0
cloud = StructArray(Obsea.State(1, 15.0, x, y, NaN, NaN) for y in xgrid, x in xgrid)
anim = @animate for t = 1:Nt
    print(t, "/", Nt, "\r")
    ℓt =
        Obsea.LikelihoodSlice(ℓ.r[:, t, :], ℓ.a[:, :, t, :], ℓ.m[t, :], ℓ.Σm[t])
    pd = models.pd
    ℓᵣ, ℓₐ = Obsea.make(ℓt, grid)
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
