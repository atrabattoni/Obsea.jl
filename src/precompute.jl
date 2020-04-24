function precompute(z, tdoalut, models, grid)
    @unpack nmode, v, τ = tdoalut
    @unpack rrange, mrange, τrange = grid
    ℓ = ones(length(rrange), length(mrange))
    for (m, model) in zip(mrange, models)
        @unpack lam = model
        u = exp.(-lam ./ 2.0) .* besseli.(0, sqrt.(lam .* z))
        for mode = 1:nmode
            A = convsame(u, v[mode])
            itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
            itp = extrapolate(itp, 1.0)
            itp = scale(itp, τrange)
            ℓ[:, m] .*= itp.(τ[:, mode])
        end
    end
    itp = interpolate(ℓ, (BSpline(Cubic(Line(OnGrid()))), NoInterp()))
    itp = extrapolate(itp, (1.0, Throw()))
    itp = scale(itp, rrange, mrange)
    (ℓ, itp)
end

function precompute(z, models, grid)
    @unpack frange, arange, mrange = grid
    ℓ = zeros(length(frange), length(arange), length(mrange))
    for (k, model) in zip(mrange, models)
        @unpack mrl, n = model
        for (j, a) in enumerate(arange)
            for i = 1:length(frange)
                ℓ[i, j, k] = wrapcauchy(z[i], a, mrl)
            end
        end
        ℓ[:, :, k] = rollprod(ℓ[:, :, k], n)
    end
    itp = interpolate(
        ℓ,
        (
            BSpline(Cubic(Line(OnGrid()))),
            BSpline(Cubic(Periodic(OnGrid()))),
            NoInterp(),
        ),
    )
    itp = extrapolate(itp, (1.0, Periodic(), Throw()))
    itp = scale(itp, frange, arange, mrange)
    (ℓ, itp)
end


function precompute(zr, za, tdoalut, models, grid)
    @unpack rrange, frange, arange, mrange = grid
    ℓr, ritp = precompute(zr, tdoalut, models, grid)
    ℓa, aitp = precompute(za, models, grid)
    itp(r, f, a, m) = ritp(r, m) * aitp(f, a, m)
    (r = ℓr, a = ℓa), itp
end
