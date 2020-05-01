import Obsea:
    TDOALUT,
    parameters,
    likelihood,
    marginalize,
    init,
    fetch,
    predict!,
    update!,
    resample!,
    estimate

using NPZ

## Load data & parameters
fs = 50.0
Np = 100000
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

@time ℓr = likelihood(zr, tdoalut, models, grid)
@time ℓa = likelihood(za, models, grid)
@time ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)

## Track
weights, particles = init(Nt, Np)
@time for t = 1:Nt
    print(t, "/", Nt, "\r")
    cloud, prevcloud, ℓ = fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    predict!(weights, cloud, prevcloud, ℓ, models, grid)
    update!(weights, cloud, ℓ, models, grid)
    resample!(weights, particles)
end
xe, ye = estimate(particles)
## Track
gr()
weights, particles = init(Nt, Np)
anim = @animate for t = 1:Nt
    print(t, "/", Nt, "\r")
    cloud, prevcloud, ℓ = fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    predict!(weights, cloud, prevcloud, ℓ, models, grid)
    update!(weights, cloud, ℓ, models, grid)
    resample!(weights, particles)
    xe, ye = estimate(particles)
    scatter(
        xe ./ 1000,
        ye ./ 1000;
        label = "estimation",
        color = "blue",
        xlim = (-30, 30),
        ylim = (-30, 30),
    )
    scatter!([x[t]] ./ 1000, [y] ./ 1000; color = "green")
end
gif(anim, "anim.gif", fps = 5)
##
using BenchmarkTools
using Profile

## Precompute likelihood
using Plots
include("../plot.jl")

## Plot true and estimated track
plotly()
scatter(
    x ./ 1000,
    fill(y / 1000, Nt);
    label = "trajectory",
    color = "green",
    xlim = (-30, 30),
    ylim = (-30, 30),
    aspect_ratio = :equal,
)
scatter!([0], [0]; label = "OBS", color = "orange")
scatter!(xe ./ 1000, ye ./ 1000; label = "estimation", color = "blue")

## Animation
weights, particles = init(Nt, Np)
gr()
anim = @animate for t = 1:Nt
    print(t, "/", Nt, "\r")
    cloud, prevcloud, ℓ = fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    predict!(weights, cloud, prevcloud, ℓ, models, grid)
    update!(weights, cloud, ℓ, models, grid)


    plot(particles[t, :], grid)

    resample!(weights, particles)
    xe, ye = estimate(particles)
    scatter!([0], [0]; color = "orange")
    scatter!([x[t]] ./ 1000, [y] ./ 1000; color = "green")
    scatter(xe ./ 1000, ye ./ 1000; label = "estimation", color = "blue")

end
gif(anim, "anim.gif", fps = 10)
xe, ye = estimate(particles)
scatter(xe ./ 1000, ye ./ 1000; label = "estimation", color = "blue")
## Plot raw data
heatmap(collect(1:Nt), grid.τ, zr; color = :viridis)
heatmap(collect(1:Nt), grid.f, za; color = :phase)

scatter(
    x ./ 1000,
    fill(y / 1000, Nt),
    color = "green",
    xlim = (-30, 30),
    ylim = (-30, 30),
    aspect_ratio = :equal,
)
scatter!([0], [0], color = "orange")

## Plot likelihood
heatmap(1:Nt, grid.r / 1000, log.(ℓr[:, :, 1]))
heatmap(1:Nt, rad2deg.(grid.a), log.(sum(ℓa, dims = 1)[1, :, :, 1]))

## Update
xgrid = -30000.0:500.0:30000.0
cloud = StructArray(State(1, 15.0, x, y, NaN, NaN) for y in xgrid, x in xgrid)
anim = @animate for t = 1:100 #
    print(t, "/", Nt, "\r")
    _, _, ℓ = fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    z = ones(size(cloud))
    update!(z, cloud, ℓ, models, grid)
    heatmap(
        xgrid / 1000,
        xgrid / 1000,
        log.(abs.(z));
        xlim = (-30, 30),
        ylim = (-30, 30),
        clim = (-20, 20),
        title = string(t),
        aspect_ratio = :equal,
    )
end
gif(anim, "anim.gif", fps = 5)

## Step by Step
weights, particles = init(Nt, Np)
t = 30
##
t += 1
cloud, prevcloud, ℓ = fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
plot(particles[t-1, :], grid)
predict!(weights, cloud, prevcloud, ℓ, models, grid)
plot!(particles[t, :], grid)
update!(weights, cloud, ℓ, models, grid)
plot!(particles[t, :], grid)
resample!(weights, particles)
plot!(particles[t, :], grid)

##
plot(ℓ.r)
heatmap(grid.a, grid.f, log.(ℓ.a[:, :, 1]))

xgrid = -30000.0:1000.0:30000.0
cloud = StructArray(State(1, 15.0, x, y, NaN, NaN) for x in xgrid, y in xgrid)
z = ones(size(cloud))
z = update!(z, cloud, ℓs[t], models, grid)
heatmap(xgrid, xgrid, log.(abs.(z)), clim = (-20, 20))

##
idx = rand(1:Np, 1000)
replace!(particles.x, 0.0 => NaN)
replace!(particles.y, 0.0 => NaN)
plot(
    particles.x[:, idx] ./ 1000,
    particles.y[:, idx] ./ 1000;
    legend = false,
    color = :black,
    alpha = 0.1,
    xlim = (-20, 20),
    ylim = (-20, 20),
    aspect_ratio = :equal,
)

## Track
gr(size = (800, 800))
weights, particles = init(Nt, Np)
anim = @animate for t = 1:Nt
    print(t, "/", Nt, "\r")
    cloud, prevcloud, ℓ = fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    predict!(weights, cloud, prevcloud, ℓ, models, grid)
    update!(weights, cloud, ℓ, models, grid)

    idx = rand(1:Np, 250)
    replace!(particles.x, 0.0 => NaN)
    replace!(particles.y, 0.0 => NaN)
    plot(
        particles.x[:, idx] ./ 1000,
        particles.y[:, idx] ./ 1000;
        legend = false,
        color = :black,
        alpha = 0.1,
        xlim = (-30, 30),
        ylim = (-30, 30),
        aspect_ratio = :equal,
        linewidth = 3.0,
    )
    scatter!(
        [x[t]] ./ 1000,
        [y] ./ 1000;
        markershape = :x,
        markersize = 10,
        markercolor = :blue,
    )

    resample!(weights, particles)
end
gif(anim, "anim.gif", fps = 5)
