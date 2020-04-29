function track(ceps, az, params, fs, N)

    # Parameters
    @assert size(ceps, 2) == size(az, 2)
    @assert isodd(size(ceps, 1))
    Nfft = (size(ceps, 1) - 1) * 2
    Nt = size(ceps, 2)
    models, propa, grid = parameters(params, fs, Nfft)
    tdoalut = TDOALUT(propa, grid)

    # Precompute

    ℓs = []
    for k = 1:Nt
        print(k, "/", Nt, "\r")
        ℓ = likelihood(ceps[:, k], az[:, k], tdoalut, models, grid)
        push!(ℓs, ℓ)
    end

    # Initialize
    weights, cloud = init(N)

    # Track
    for k = 1:Nt
        predict!(weights, cloud, ℓs[k], models, grid)
        update!(weights, cloud, ℓs[k], models, grid)
        weights, cloud = resample(weights, cloud)
    end

    # Finalize
    final!(cloud)

    # Output
    @unpack Nm = grid
    estimate(cloud, Nm, Nt)
end
