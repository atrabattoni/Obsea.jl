function track(zr, za, params, fs, Nfft, Np)

    # Parameters
    @assert size(zr, 2) == size(za, 2)
    Nt = size(zr, 2)
    models, propa, grid = parameters(params, fs, Nfft)
    tdoalut = TDOALUT(propa, grid)

    # Precompute
    ℓr = likelihood(zr, tdoalut, models, grid)
    ℓa = likelihood(za, models, grid)
    ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)

    # Track
    weights, particles = init(Nt, Np)
    for t = 1:Nt
        cloud, prevcloud, ℓ = fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
        predict!(weights, cloud, prevcloud, ℓ, models, grid)
        update!(weights, cloud, ℓ, models, grid)
        resample!(weights, particles)
    end

    # Output
    return estimate(particles)
end

function fetch(t, particles, ℓr, ℓa, ℓm, ℓΣm)
    kprev = (t == 1 ? 1 : t - 1)
    @views ℓ = Likelihood(ℓr[:, t, :], ℓa[:, :, t, :], ℓm[t, :], ℓΣm[t])
    @views cloud = particles[t, :]
    @views prevcloud = particles[kprev, :]
    return cloud, prevcloud, ℓ
end
