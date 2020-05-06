"""
    Model(name, q, vmin, vmax, ps, pb, pd, lam, mrl, n)

Store all parameters related to an object model.

> Note: all parameters must be provided as keyword arguments.

# Arguments
- `name::String`: the name of the model.
- `q::Float64`: the process noise (m/s²).
- `vmin::Float64`: the minimum speed (m/s).
- `vmax::Float64`: the maximum speed (m/s).
- `ps::Float64`: the survival probability.
- `pb::Float64`: the birth probability.
- `pd::Float64`: the detection probability.
- `lam::Float64`: the cepstral SNR.
- `mrl::Float64`: the mean running length.
- `n::Int`: the sectral width (px).
"""
@with_kw struct Model
    name::String
    q::Float64
    vmin::Float64
    vmax::Float64
    ps::Float64
    pb::Float64
    pd::Float64
    lam::Float64
    mrl::Float64
    n::Int
end

"""
    Propagation(Nmode, depth, celerity, ic, ib, sigma)

Store all parameters related to wave propagation.

> Note: all parameters must be provided as keyword arguments.

# Arguments
- `Nmode::Int`: the number of interfering paths.
- `depth::Float64`: the ocean depth (m).
- `celerity::Foat64`: the sound equivalent speed in the Ocean (m/s).
- `ic::Float64`: the critical incidence angle (°).
- `ib::Float64`: the bouncing incidence angle (°).
- `sigma::Vector{Float64}` the quefrency incertitude for each mode (px).
"""
struct Propagation
    Nmode::Int
    depth::Float64
    celerity::Float64
    ic::Float64
    ib::Float64
    sigma::Vector{Float64}
    function Propagation(; Nmode, depth, celerity, ic, ib, sigma)
        ic = deg2rad(ic)
        ib = deg2rad(ib)
        new(Nmode, depth, celerity, ic, ib, sigma)
    end
end

"""
    Grid(Nm, fs, nfft; rmax, rres, fmin, fmax, ares)

Store all parameters related to the computational grid.

# Arguments
- `Nm::Int`: the number of models.
- `fs::Float64`: the sampling rate (Hz).
- `nfft::Int`: the FFT window size.
- `rmax::Float64`: the maximum range (m).
- `rres::Float64`: the range spacing (m).
- `fmin::Float64`: the minimum frequency (Hz).
- `fmax::Float64`: the maximum frequency (Hz).
- `ares::Int`: the number of azimuth grid points.
"""
struct Grid
    r::Vector{Float64}
    τ::Vector{Float64}
    f::Vector{Float64}
    a::Vector{Float64}
    Nr::Int
    Nτ::Int
    Nf::Int
    Na::Int
    Nm::Int
    T::Float64
    function Grid(Nm, fs, nfft; rmax, rres, fmin, fmax, ares)
        r = collect(0.0:rres:rmax)
        τ = collect(range(0, nfft / fs / 2, length = (nfft ÷ 2) + 1))
        f = collect(range(0, fs / 2, length = (nfft ÷ 2) + 1))
        f = limit(f, fmin, fmax)
        a = collect(range(0, 2π, length = ares + 1))
        Nr = length(r)
        Nτ = length(τ)
        Nf = length(f)
        Na = length(a)
        T = nfft / fs / 2
        new(r, τ, f, a, Nr, Nτ, Nf, Na, Nm, T)
    end
end

"""
    parameters(params, fs, nfft)

Parse parameters from a dictionnary or a TOML file.
"""
function parameters(dict::Dict, fs, nfft)
    models = StructVector([Model(; symbolize(d)...) for d in dict["model"]])
    propa = Propagation(; symbolize(dict["propagation"])...)
    grid = Grid(length(models), fs, nfft; symbolize(dict["grid"])...)
    return (models, propa, grid)
end
function parameters(fname::String, fs, nfft)
    dict = TOML.parsefile(fname)
    return parameters(dict::Dict, fs, nfft)
end
