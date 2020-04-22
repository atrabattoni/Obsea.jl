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


@with_kw struct Propagation
    nmode
    depth
    celerity
    ic
    ib
    sigma
end


struct Grid
    τrange
    rrange
    frange
    arange
    mrange
    T
    function Grid(nmodel, fs, nfft; rmax, rres, fmin, fmax, ares)
        τrange = range(0, nfft / fs / 2, length = nfft ÷ 2 + 1)
        rrange = range(0, rmax, step = rres)
        frange = range(0, fs / 2, length = nfft ÷ 2 + 1)
        frange = limit(frange, fmin, fmax)
        arange = range(0, 2π, length = ares + 1)
        mrange = range(1, nmodel, step = 1)
        T = nfft / fs / 2
        new(τrange, rrange, frange, arange, mrange, T)
    end
end


function parameters(dict::Dict, fs, nfft)
    models = [Model(; symbolize(d)...) for d in dict["model"]]
    nmodel = propa = Propagation(; symbolize(dict["propagation"])...)
    grid = Grid(length(models), fs, nfft; symbolize(dict["grid"])...)
    return (models, propa, grid)
end

function parameters(fname::String, fs, nfft)
    dict = TOML.parsefile(fname)
    return parameters(dict::Dict, fs, nfft)
end
