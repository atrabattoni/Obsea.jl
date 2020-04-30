function track(ceps, az, params, fs, Nfft, Np)

    # Parameters
    @assert size(ceps, 2) == size(az, 2)
    Nt = size(ceps, 2)
    models, propa, grid = parameters(params, fs, Nfft)
    tdoalut = TDOALUT(propa, grid)

    # Precompute
    @views ℓs =
        [likelihood(ceps[:, k], az[:, k], tdoalut, models, grid) for k = 1:Nt]

    # Track
    weights, particles = init(Nt, Np)
    for k = 1:Nt
        cloud, prevcloud, ℓ = fetch(k, particles, ℓs)
        predict!(weights, cloud, prevcloud, ℓ, models, grid)
        update!(weights, cloud, ℓ, models, grid)
        resample!(weights, particles)
    end

    # Output
    return estimate(particles)
end

function fetch(k, particles, ℓs)
    kprev = (k == 1 ? 1 : k - 1)
    ℓ = ℓs[k]
    @views cloud = particles[k, :]
    @views prevcloud = particles[kprev, :]
    return cloud, prevcloud, ℓ
end
