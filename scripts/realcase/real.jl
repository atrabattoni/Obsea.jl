import Obsea:
    TDOALUT,
    parameters,
    likelihood,
    rollprod,
    init,
    predict!,
    update!,
    resample,
    final!,
    estimate,
    marginalize,
    interp,
    State,
    grid

using NPZ
using Plots
using BenchmarkTools
using Profile
include("plot.jl")

## Load data & parameters
fs = 50.0
Np = 5000

ceps = npzread("ceps.npy")
az = npzread("az.npy")

Nfft = (size(ceps, 1) - 1) * 2
Nt = size(ceps, 2)

models, propa, grid = parameters("params.toml", fs, Nfft)
tdoalut = TDOALUT(propa, grid)

## Plot data
heatmap(collect(1:Nt), grid.τ, ceps; color = :viridis)
heatmap(collect(1:Nt), grid.f, az; color = :phase)

## Precompute likelihood
ℓs = []
for k = 1:Nt
    print(k, "/", Nt, "\r")
    ℓ = likelihood(ceps[:, k], az[:, k], tdoalut, models, grid)
    push!(ℓs, ℓ)
end

ℓr = hcat([ℓ.r for ℓ in ℓs]...)
heatmap(1:Nt, grid.r / 1000, log.(ℓr))

ℓa = hcat([sum(ℓ.a, dims = (1, 3))[1, :, 1] for ℓ in ℓs]...)
heatmap(1:Nt, rad2deg.(grid.a), log.(ℓa))

## Track

function iterate!(weights, cloud, ℓs, models, grid)
    for k = 1:100 # anim = @animate
        print(k, "/", Nt, "\r")
        ℓ =
            ℓs[k]::NamedTuple{
                (:r, :a, :m, :Σm),
                Tuple{
                    Array{Float64,2},
                    Array{Float64,3},
                    Array{Float64,1},
                    Float64,
                },
            }
        predict!(weights, cloud, ℓ, models, grid)
        update!(weights, cloud, ℓ, models, grid)
        weights, cloud = resample(weights, cloud)
        # plot(cloud, grid)
    end
    # anim
end

weights, cloud = init(Np)
@code_warntype iterate!(weights, cloud, ℓs, models, grid)

@profiler iterate!(weights, cloud, ℓs, models, grid)
@time iterate!(weights, cloud, ℓs, models, grid)
gif(anim, "anim.gif", fps = 10)
##
@profiler weights, cloud = resample(weights, cloud)

# Finalize
final!(cloud)

# Output
m, x, y = estimate(cloud, grid.Nm, Nt)

plot(m')
plot(x, y, xlim = (-30000, 30000), ylim = (-30000, 30000))
## Update

anim = @animate for k = 1:100
    print(k, "/", Nt, "\r")
    xgrid = -30000.0:1000.0:30000.0
    z = [
        interp(ℓs[k], State(1, 15.0, x, y, NaN, NaN), models, grid)
        for x in xgrid, y in xgrid
    ]
    heatmap(xgrid, xgrid, log.(abs.(z)), clim = (-20, 20))
end
gif(anim, "anim.gif", fps = 10)
