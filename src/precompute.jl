function precompute(z, tdoalut, models, grid)
    @unpack Nmode, v, τ = tdoalut
    @unpack Nr, Nm = grid
    ℓ = ones(Nr, Nm)
    for m = 1:Nm
        @unpack lam = models[m]
        u = exp.(-lam ./ 2.0) .* besseli.(0, sqrt.(lam .* z))
        for mode = 1:Nmode
            A = convsame(u, v[mode])
            itp = interpolate(A, BSpline(Cubic(Line(OnGrid())))) #
            itp = extrapolate(itp, 1.0)
            itp = scale(itp, grid2range(grid.τ))
            ℓ[:, m] .*= itp.(τ[:, mode])
        end
    end

    ℓ
end

function precompute(z, models, grid)
    @unpack Nf, Na, Nm = grid
    ℓ = zeros(Nf, Na, Nm)
    for k = 1:Nm
        @unpack mrl, n = models[k]
        for (j, a) in enumerate(grid.a)
            for i = 1:Nf
                ℓ[i, j, k] = wrapcauchy(z[i], a, mrl)
            end
            ℓ[:, j, k] = rollprod(ℓ[:, j, k], n)
        end
    end
    ℓ
end


function precompute(zr, za, tdoalut, models, grid)
    ℓr = precompute(zr, tdoalut, models, grid)
    ℓa = precompute(za, models, grid)
    (r = ℓr, a = ℓa)
end
