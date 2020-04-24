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


function birth(ℓ, models, grid)
    @unpack rrange, frange, arange, mrange = grid
    N = length(rrange) * length(frange) * length(arange)
    pbs = [model.pb for model in models]
    normalization =
        (1.0 - sum(pbs)) + sum(pbs[m] * sum(ℓ[:, :, :, m]) / N for m in mrange)
    p = cat([pbs[m] * ℓ[:, :, :, m] / N for m in mrange]..., dims = 4)
    cdf = cumsum(vec(p)) / normalization
    idx = argsample(cdf)
    if !iszero(idx)
        idx = CartesianIndices(p)[idx]
        r = rrange[idx[1]]
        f = frange[idx[2]]
        a = arange[idx[3]]
        m = mrange[idx[4]]
        @unpack vmin, vmax = models[m]
        vr = vmin + (vmax - vmin) * rand()
        va = 2π * rand()
        x = r * sin(a)
        y = r * cos(a)
        vx = vr * sin(va)
        vy = vr * cos(va)
        return (normalization / ℓ[idx], State(m, f, x, y, vx, vy))
    else
        return (normalization, State())
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
