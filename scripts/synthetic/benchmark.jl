using NPZ
import Obsea:
    parameters,
    TDOALUT,
    likelihood,
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

zr = npzread("ceps_test.npy")
za = npzread("az_test.npy")

fs = 50.0
Np = 10000
Nfft = 1024
models, propa, grid = parameters("params_test.toml", 50.0, 1024)


function precompute(zr, za, models, propa, grid)
    tdoalut = TDOALUT(propa, grid)
    ℓr = likelihood(zr, tdoalut, models, grid)
    ℓa = likelihood(za, models, grid)
    ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)
    return ℓr, ℓa, ℓm, ℓΣm
end

println("precompute")
ℓr, ℓa, ℓm, ℓΣm = precompute(zr, za, models, propa, grid)
@btime precompute(zr, za, models, propa, grid)
@profiler precompute(zr, za, models, propa, grid)

weights, particles = init(Nt, Np)

function loop(weights, particles, ℓr, ℓa, ℓm, ℓΣm, models, grid, Np)
    Nt = size(ℓr, 2)

    for t = 1:Nt
        cloud, prevcloud, ℓ = slice(t, particles, ℓr, ℓa, ℓm, ℓΣm)
        predict!(weights, cloud, prevcloud, ℓ, models, grid)
        update!(weights, cloud, ℓ, models, grid)
        if (1.0 / sum(weights .^ 2)) < (Np / 2.0)
            resample!(weights, particles)
        end
    end
end

println("loop")
weights, particles = init(Nt, Np)
loop(weights, particles, ℓr, ℓa, ℓm, ℓΣm, models, grid, Np)
@btime loop(w, p, ℓr, ℓa, ℓm, ℓΣm, models, grid, Np) setup =
    ((w, p) = init(Nt, Np)) evals = 1
weights, particles = init(Nt, Np)
@profiler loop(weights, particles, ℓr, ℓa, ℓm, ℓΣm, models, grid, Np)

Nt = size(ℓr, 2)

particles = StructArray{State}(undef, Nt, Np)
@btime fill!(particles, State())
