struct CDF
    r::Array{Float64,1}
    a::Array{Float64,1}
end


struct ITP
    r::ScaledInterpolation
    a::ScaledInterpolation
end


function precompute(z, tdoalut, models, grid)
    @unpack nmode, v, τ = tdoalut
    @unpack rrange, mrange, τrange = grid
    y = ones(length(rrange), length(mrange))
    for (m, model) in zip(mrange, models)
        @unpack lam = model
        u = exp.(-lam ./ 2.0) .* besseli.(0, sqrt.(lam .* z))
        for mode = 1:nmode
            A = convsame(u, v[mode])
            itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
            itp = extrapolate(itp, 1.0)
            itp = scale(itp, τrange)
            y[:, m] .*= itp.(τ[:, mode])
        end
    end

    cdf = cumsum(vec(y))
    cdf /= last(cdf)

    itp = interpolate(y, (BSpline(Cubic(Line(OnGrid()))), NoInterp()))
    itp = extrapolate(itp, (1.0, Throw()))
    itp = scale(itp, rrange, mrange)

    (cdf, itp)
end

function precompute(z, models, grid)
    @unpack frange, arange, mrange = grid
    y = zeros(length(frange), length(arange), length(mrange))
    for (k, model) in zip(mrange, models)
        @unpack mrl, n = model
        for (j, a) in enumerate(arange)
            for i = 1:length(frange)
                y[i, j, k] = wrapcauchy(z[i], a, mrl)
            end
        end
        y[:, :, k] = rollprod(y[:, :, k], n)
    end

    cdf = cumsum(vec(y))
    cdf /= last(cdf)

    itp = interpolate(
        y,
        (
            BSpline(Cubic(Line(OnGrid()))),
            BSpline(Cubic(Periodic(OnGrid()))),
            NoInterp(),
        ),
    )
    itp = extrapolate(itp, (1.0, Periodic(), Throw()))
    itp = scale(itp, frange, arange, mrange)

    (cdf, itp)
end

function precompute(zr, za, tdoalut, models, grid)
    rcdf, ritp = precompute(zr, tdoalut, models, grid)
    acdf, aitp = precompute(za, models, grid)
    (CDF(rcdf, acdf), ITP(ritp, aitp))
end
