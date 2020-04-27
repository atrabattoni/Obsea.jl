function predict!(weights, cloud, ℓ, models, grid)
    @assert length(weights) == length(cloud)
    for (i, particle) in enumerate(cloud)
        normalization, state = transition(last(particle), ℓ, models, grid)
        weights[i] *= normalization
        push!(particle, state)
    end
end


function transition(state, ℓ, models, grid)
    if isempty(state)
        return birth(ℓ, models, grid)
    else
        return move(state, models, grid)
    end
end


function move(state, models, grid)
    @unpack T = grid
    @unpack ps, q = models[getmodel(state)]
    @unpack m, f, x, y, vx, vy = state
    if rand() < ps
        ax = q * randn()
        ay = q * randn()
        x = x + vx * T + ax * T^2 / 2.0
        y = y + vy * T + ay * T^2 / 2.0
        vx = vx + ax * T
        vy = vy + ay * T
        return (1.0, State(m, f, x, y, vx, vy))
    else
        return (1.0, State())
    end
end


function randspeed(model)
    @unpack vmin, vmax = model
    vr = vmin + (vmax - vmin) * rand()
    va = 2π * rand()
    va, vr
end


function birth(ℓ, models, grid)
    pb = [model.pb for model in models]
    @unpack Nr, Na, Nf, Nm, rrange, frange, arange, mrange = grid
    ℓ0 = 1.0 - sum(pb)
    ℓm =
        [pb[m] * sum(ℓ.r[:, m]) / Nr * sum(ℓ.a[:, :, m]) / Na / Nf for m = 1:Nm]
    normalization = ℓ0 + sum(ℓm)
    # model
    cdf = cumsum(ℓm)
    m = argsample(cdf, scale = normalization)
    if !iszero(m)
        # range
        cdf = cumsum(ℓ.r[:, m])
        @show first(cdf) last(cdf) last(cdf) + 10 * eps(last(cdf))
        idx = argsample(cdf, scale = last(cdf))
        normalization /= ℓ.r[idx, m]
        r = rrange[idx]
        # frequency & azimuth
        cdf = cumsum(vec(ℓ.a[:, :, m]))
        idx = argsample(cdf, scale = last(cdf))
        idx = CartesianIndices(ℓ.a[:, :, m])[idx]
        normalization /= ℓ.a[idx[1], idx[2], m]
        f = frange[idx[1]]
        a = arange[idx[2]]
        vr, va = randspeed(models[m])
        x, y, = ra2xy(r, a)
        vx, vy = ra2xy(vr, va)
        return normalization, State(m, f, x, y, vx, vy)
    else
        return normalization, State()
    end
end


function logf(state, prevstate, models, grid)
    @unpack T = grid
    if isempty(state)
        if isempty(prevstate)
            return log(1.0 - pb)
        else
            return log(1.0 - ps)
        end
    else
        if isempty(prevstate)
            return log(pb)
        else
            @assert getmodel(state) == getmodel(prevstate)
            @unpack q, pb, ps = models[getmodel(state)]
            dvx = state.vx - prevstate.vx
            dvy = state.vy - prevstate.vy
            return log(ps) - (dvx^2 + dvy^2) / (q * T)^2  # TODO
        end
    end
end
