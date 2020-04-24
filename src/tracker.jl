function track(ceps, az, params, fs, N)

    # Parameters
    @assert size(ceps) == size(az)
    @assert isodd(size(ceps, 1))
    nfft = (size(ceps, 1) - 1) * 2
    nt = size(ceps, 2)
    models, propa, grid = parameters(params, fs, nfft)

    # Precompute
    ys = Array{Float64,5}(
        undef,
        length(rrange),
        length(frange),
        length(arange),
        length(mrangge),
        nt,
    )
    itps = Vector{ScaledInterpolation}(undef, nt)
    for k = 1:nt
        ys[k], itps[k] = precompute(ceps[:, k], az[:, k], propa, models, grid)
    end

    # Init
    weights, cloud = init(N)

    # Track
    for k = 1:nt
        predict!(weights, cloud, cdfs[k], models, grid)
        update!(weights, cloud, itps[k], models)
        resample!(weights, cloud)
    end

    # Output
    cloud
end
