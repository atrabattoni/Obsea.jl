function precompute(z, tdoalut, model, grid)
    @unpack nmode, v, τ = tdoalut
    @unpack lam = model
    @unpack τrange = grid
    u = exp.(-lam ./ 2.0) .* besseli.(0, sqrt.(lam .* z))
    y = ones(size(τ, 1))
    for mode = 1:nmode
        A = convsame(u, v[mode])
        itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
        itp = extrapolate(itp, 1.0)
        itp = scale(itp, τrange)
        y .*= itp(τ[:, mode])
    end

    cdf = cumsum(y)
    cdf /= last(cdf)

    itp = interpolate(y, BSpline(Cubic(Line(OnGrid()))))
    itp = extrapolate(itp, 1.0)
    itp = scale(itp, rrange)

    (cdf, itp)
end


function precompute(z, models, grid)
    @unpack frange, arange, mrange = grid
    y = zeros(length(frange), length(arange), length(mrange))
    for (k, model) in zip(mrange, models)
        @unpack mrl = model
        for j = 1:length(arange)
            for i = 1:length(frange)
                y[i, j, k] = wrapcauchy(z[i], a[j], mrl)
            end
        end
        y[:, :, k] = rollprod([:, :, k], n)
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


struct Cdf
    r::Array{Float64,1}
    fam::Array{Float64,1}
end


struct Itp
    r::ScaledInterpolation
    fam::ScaledInterpolation
end


function precompute(zr, zfam, tdoalut, models, grid)
    rcdf, citp = precompute(zr, tdoalut, models, grid)
    famcdf, famitp = precompute(zfam, models, grid)
    (Cdf(rcdf, famcdf), Itp(ritp, famitp))
end
