@with_kw struct Model
    name
    q
    vmin
    vmax
    ps
    pb
    pd
    lam
    mrl
    n
end


struct Propagation
    Nmode
    depth
    celerity
    ic
    ib
    sigma
    function Propagation(; Nmode, depth, celerity, ic, ib, sigma)
        ic = deg2rad(ic)
        ib = deg2rad(ib)
        new(Nmode, depth, celerity, ic, ib, sigma)
    end
end


struct Grid
    τrange
    rrange
    frange
    arange
    mrange
    Nτ
    Nr
    Nf
    Na
    Nm
    T
    function Grid(Nm, fs, nfft; rmax, rres, fmin, fmax, ares)
        τrange = range(0, nfft / fs / 2, length = nfft ÷ 2 + 1)
        rrange = range(0, rmax, step = rres)
        frange = range(0, fs / 2, length = nfft ÷ 2 + 1)
        frange = limit(frange, fmin, fmax)
        arange = range(0, 2π, length = ares + 1)
        mrange = range(1, Nm, step = 1)
        Nτ = length(τrange)
        Nr = length(rrange)
        Nf = length(frange)
        Na = length(arange)
        Nm = length(mrange)
        T = nfft / fs / 2
        new(τrange, rrange, frange, arange, mrange, Nτ, Nr, Nf, Na, Nm, T)
    end
end


function parameters(dict::Dict, fs, nfft)
    models = [Model(; symbolize(d)...) for d in dict["model"]]
    Nm = propa = Propagation(; symbolize(dict["propagation"])...)
    grid = Grid(length(models), fs, nfft; symbolize(dict["grid"])...)
    return (models, propa, grid)
end

function parameters(fname::String, fs, nfft)
    dict = TOML.parsefile(fname)
    return parameters(dict::Dict, fs, nfft)
end
