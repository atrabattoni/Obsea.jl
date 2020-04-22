import pkg.TOML.parsefile


struct Propagation
    nmode::Int64
    depth::Float64
    celerity::Float64
    ic::Float64
    ib::Float64
    lam::Float64
    v::Dict{Int64,Array{Float64,1}}
    τ::Array{Float64,2}
    function Propagation(; nmode, depth, celerity, ic, ib, lam, n, σ)
        @assert length(n) == length(σ) == nmode
        @assert depth > 0
        @assert celerity > 0
        @assert 0 <= ic
        v = Dict(mode => window(n[mode], σ[mode]) for mode = 1:length(n))
        τ = [
            tdoa(r, mode, depth, celerity, ic, ib)
            for r in rrange, mode = 1:nmode
        ]
        new(nmode, depth, celerity, ic, ib, lam, v)
    end
end


@with_kw struct Model
    q::Float64 = 0.0
    vmin::Float64 = 0.0
    vmax::Float64 = Inf
    ps::Float64 = 1.0
    pb::Float64 = 0.0
    pd::Float64 = 1.0
    n::Int64 = 1
end


struct Grid
    τrange::AbstractRange
    rrange::AbstractRange
    frange::AbstractRange
    arange::AbstractRange
    mrange::UnitRange{Int64}
    T::Float64
    function Grid(d, nmodel, fs, nfft)
        τrange = range(0, nfft / fs / 2, length = nfft ÷ 2 + 1)
        rrange = range(0, d["rmax"], step = d["rres"])
        frange = range(0, fs / 2, length = nfft ÷ 2 + 1)
        frange = limit(frange, d["fmin"], d["fmax"])
        arange = range(0, 2π, length = d["na"] + 1)
        mrange = range(1, nmodel)
        T = nfft / fs / 2
        new(τrange, rrange, frange, arange, mrange, T)
    end
end


function parameters(dict::Dict, fs, nfft)
    models = [Model(; symbolize(d)...) for d in dict["model"]]
    propa = Propagation(; symbolize(dict["range"])...)
    grid = Grid(dict["grid"], length(dict["model"]), fs, nfft)
    return (models, propa, grid)
end

function parameters(fname::String, fs, nfft)
    dict = parsefile(fname)
    return parameters(dict::Dict, fs, nfft)
end
