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
    for m = 1:Nm
        @unpack mrl, n = models[m]
        for (j, a) in enumerate(grid.a)
            for i = 1:Nf
                ℓ[i, j, m] = wrapcauchy(z[i], a, mrl)
            end
            ℓ[:, j, m] = rollprod(ℓ[:, j, m], n)
        end

    end
    ℓ
end

function precompute(zr, za, tdoalut, models, grid)
    ℓr = precompute(zr, tdoalut, models, grid)
    ℓa = precompute(za, models, grid)
    ℓ = (r = ℓr, a = ℓa)
end


function distribution(ℓ, models, grid)
    @unpack Nr, Nf, Na, Nm = grid
    # model
    pb = [model.pb for model in models]
    ℓm = similar(pb)
    for m = 1:Nm
        @views ℓm[m] = pb[m] * sum(ℓ.r[:, m]) / Nr * sum(ℓ.a[:, :, m]) / Na / Nf
    end
    Fm = cumsum(ℓm)
    ℓ0 = 1.0 - sum(pb)
    FΣm = ℓ0 + last(Fm)
    # range
    Fr = similar(ℓ.r)
    for m = 1:Nm
        @views Fr[:, m] = cumsum(ℓ.r[:, m])
    end
    # azimuth
    Fa = Array{Float64,2}(undef, Nf * Na, Nm)
    for m = 1:Nm
        @views Fa[:, m] = cumsum(vec(ℓ.a[:, :, m]))
    end
    F = (r = Fr, a = Fa, m = Fm, Σm = FΣm)
end
