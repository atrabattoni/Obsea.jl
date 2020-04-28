function predict!(weights, cloud, ℓ, F, models, grid)
    mask = isempty(cloud)
    birth!(view(weights, mask), view(cloud, mask), ℓ, F, models, grid)
    move!(view(cloud, .!mask), models, grid)
end

# for (i, particle) in enumerate(cloud)
#     normalization, state = transition(last(particle), ℓ, F, models, grid)
#     weights[i] *= normalization
#     push!(particle, state)
# end

# function transition(state, ℓ, F, models, grid)
#     if isempty(state)
#         return birth(ℓ, F, models, grid)
#     else
#         return move(state, models, grid)
#     end
# end

function move!(cloud, models, grid)
    for particle in cloud
        state = move(last(particle), models, grid)
        push!(particle, state)
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
    else
        return State()
    end
end


function randspeed(model)
    @unpack vmin, vmax = model
    vr = vmin + (vmax - vmin) * rand()
    va = 2π * rand()
    va, vr
end


function birth!(weights, cloud, ℓ, F, models, grid)
    for (i, particle) in enumerate(cloud)
        normalization, state = birth(ℓ, F, models, grid)
        weights[i] *= normalization
        push!(particle, state)
    end
end

function birth(ℓ, F, models, grid)
    @unpack Nr, Na, Nf, Nm = grid
    # model
    normalization = F.Σm
    m = argsample(F.m, scale = F.Σm)
    if !iszero(m)
        # range
        idx = argsample(F.r[:, m], scale = last(F.r[:, m]))
        normalization /= ℓ.r[idx, m]
        r = grid.r[idx]
        # frequency & azimuth
        idx = argsample(F.a[:, m], scale = last(F.a[:, m]))
        idx = CartesianIndices(ℓ.a[:, :, m])[idx]
        normalization /= ℓ.a[idx[1], idx[2], m]
        f = grid.f[idx[1]]
        a = grid.a[idx[2]]
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
