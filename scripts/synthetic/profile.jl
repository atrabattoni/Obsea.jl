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
using Profile

## Load data & parameters

Nτ = 513
Nf = 266
Nt = 100
Na = 361
Nm = 2

zr = rand(Nτ, Nt)
za = 2π .* rand(Nf, Nt)
a = collect(range(0, 2π, length = Na))
mrl = rand(Nm)

fs = 50.0
Np = 10
Nfft = 1024
models, propa, grid = parameters("params_test.toml", fs, Nfft)





## Compile
xe, ye = track(zr, za, models, propa, grid, Np)

## Check type stability
@code_warntype track(zr, za, models, propa, grid, Np)

## Benchmark
@time xe, ye = track(zr, za, models, propa, grid, Np)

## Profile
@profiler xe, ye = track(zr, za, models, propa, grid, Np)
