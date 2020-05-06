function track(zr, za, models, propa, grid, Np)

    # Parameters
    @assert size(zr, 2) == size(za, 2)
    Nt = size(zr, 2)

    # Precompute
    tdoalut = TDOALUT(propa, grid)
    ℓr = likelihood(zr, tdoalut, models, grid)
    ℓa = likelihood(za, models, grid)
    ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)
    ℓ = Likelihood(ℓr, ℓa, ℓm, ℓΣm)

    # Track
    weights, particles = init(Nt, Np)
    for t = 1:Nt
        cloud, prevcloud, ℓt = slice(t, particles, ℓ)
        predict!(weights, cloud, prevcloud, ℓt, models, grid)
        update!(weights, cloud, ℓt, models, grid)
        resample!(weights, particles, ℓ, models, grid)
    end

    # Output
    return estimate(particles)
end

function slice(t, particles, ℓ)
    kprev = (t == 1 ? 1 : t - 1)
    @views ℓt =
        LikelihoodSlice(ℓ.r[:, t, :], ℓ.a[:, :, t, :], ℓ.m[t, :], ℓ.Σm[t])
    @views cloud = particles[t, :]
    @views prevcloud = particles[kprev, :]
    return cloud, prevcloud, ℓt
end
