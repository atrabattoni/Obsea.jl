function predict!(cloud, cdf, models, grid)
    for particle in cloud
        push!(particle, transition(last(particle), cdf, models, grid))
    end
end


function transition(state, cdf, models, grid)
    if isempty(state)
        return birth(cdf, models, grid)
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
        return State(m, f, x, y, vx, vy)
    else:
        return State()
    end
end


function birth(cdf, models, grid)
    @unpack rrange, frange, arange, mrange = grid
    m = sample(cumsum(model.ps for model in models))
    if m in mrange
        idx = argsample(cdf.r)
        r = rrange[idx]

        idx = argsample(cdf.fam)
        idx = CartesianIndex((length(frange), length(arange), length(mrange)))
        f = frange[idx[1]]
        a = arange[idx[2]]
        m = mrange[idx[3]]

        @unpack vmin, vmax = models[m]
        vr = vmin + (vmax - vmin) * rand()
        va = 2π * rand()

        x = r * sin(a)
        y = r * cos(a)
        vx = vr * sin(va)
        vy = vr * cos(va)

        return State(m, f, x, y, vx, vy)
    else
        return State()
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
