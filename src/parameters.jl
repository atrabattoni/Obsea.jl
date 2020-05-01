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
