function track(zr, za, models, propa, grid, Np)

    # Parameters
    @assert size(zr, 2) == size(za, 2)
    Nt = size(zr, 2)

    # Precompute
    tdoalut = TDOALUT(propa, grid)
    ℓr = likelihood(zr, tdoalut, models, grid)
    ℓa = likelihood(za, models, grid)
    ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)

    # Track
    weights, particles = init(Nt, Np)
    for t = 1:Nt
        cloud, prevcloud, ℓt = slice(t, particles, ℓr, ℓa, ℓm, ℓΣm)
        predict!(weights, cloud, prevcloud, ℓt, models, grid)
        update!(weights, cloud, ℓt, models, grid)
        resample!(weights, particles)
    end

    # Output
    return estimate(particles)
end

function slice(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    kprev = (t == 1 ? 1 : t - 1)
    @views ℓt = Likelihood(ℓr[:, t, :], ℓa[:, :, t, :], ℓm[t, :], ℓΣm[t])
    @views cloud = particles[t, :]
    @views prevcloud = particles[kprev, :]
    return cloud, prevcloud, ℓt
end
