using NPZ
import Obsea:
    parameters,
    TDOALUT,
    likelihood,
    Likelihood,
    marginalize,
    init,
    slice,
    predict!,
    update!,
    resample!,
    predict!,
    track
using BenchmarkTools, Profile

## Load data & parameters

zr = npzread("data/ceps.npy")
za = npzread("data/az.npy")

fs = 50.0
Np = 10000
Nfft = 1024
models, propa, grid = parameters("params.toml", 50.0, 1024)
tdoalut = TDOALUT(propa, grid)
Nt = size(zr, 2)


println("precompute")
ℓ = Likelihood(zr, za, tdoalut, models, grid)
@btime Likelihood(zr, za, tdoalut, models, grid)
@profiler Likelihood(zr, za, tdoalut, models, grid)

weights, particles = init(Nt, Np)


function loop(weights, particles, ℓ, models, grid, Np, Nt)
    for t = 1:Nt
        cloud, prevcloud, ℓt = slice(t, particles, ℓ)
        predict!(weights, cloud, prevcloud, ℓt, models, grid)
        update!(weights, cloud, ℓt, models, grid)
        @views resample!(weights, particles[1:t, :], ℓ, models, grid)
    end
end

println("loop")
weights, particles = init(Nt, Np)
loop(weights, particles, ℓ, models, grid, Np, Nt)
@btime loop(w, p, ℓ, models, grid, Np, Nt) setup = ((w, p) = init(Nt, Np)) evals =
    1
weights, particles = init(Nt, Np)
@profiler loop(weights, particles, ℓ, models, grid, Np, Nt)

particles = StructArray{State}(undef, Nt, Np)
@btime fill!(particles, State())
