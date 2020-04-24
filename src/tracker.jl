function track(ceps, az, params, fs, N)

    # Parameters
    @assert size(ceps, 2) == size(az, 2)
    @assert isodd(size(ceps, 1))
    nfft = (size(ceps, 1) - 1) * 2
    nt = size(ceps, 2)
    models, propa, grid = parameters(params, fs, nfft)
    tdoalut = TDOALUT(propa, grid)

    # Precompute

    ℓs = []
    itps = []
    for k = 1:nt
        ℓ, itp = precompute(ceps[:, k], az[:, k], tdoalut, models, grid)
        push!(ℓs, ℓ)
        push!(itps, itp)
    end

    # Init
    weights, cloud = init(N)

    # Track
    for k = 1:nt
        predict!(weights, cloud, ℓs[k], models, grid)
        update!(weights, cloud, itps[k], models)
        weights, cloud = resample(weights, cloud)
    end

    # Output
    cloud
end
