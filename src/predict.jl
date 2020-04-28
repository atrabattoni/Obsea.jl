function predict!(weights, cloud, ℓ, F, models, grid)
    mask = isempty(cloud)
    birth!(view(weights, mask), view(cloud, mask), ℓ, F, models, grid)
    move!(view(cloud, .!mask), models, grid)
end


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


function counts(x, N)
    d = Dict(i => 0 for i = 0:N)
    for xi in x
        d[xi] += 1
    end
    d
end


function randspeed(model)
    @unpack vmin, vmax = model
    vr = vmin + (vmax - vmin) * rand()
    va = 2π * rand()
    va, vr
end


function birth(m, idx, ℓ, models, grid)
    value = ℓ.r[idx.r, m] * ℓ.a[idx.f, idx.a, m]
    r = grid.r[idx.r]
    f = grid.f[idx.f]
    a = grid.a[idx.a]
    vr, va = randspeed(models[m])
    x, y, = ra2xy(r, a)
    vx, vy = ra2xy(vr, va)
    return value, State(m, f, x, y, vx, vy)
end

function birth!(weights, cloud, ℓ, F, models, grid)
    @unpack Nr, Na, Nf, Nm = grid

    ms = argsample(F.m, length(cloud), scale = F.Σm)
    N = counts(ms, Nm)

    idxs = []
    @views for m = 1:Nm
        ridx = argsample(F.r[:, m], N[m], scale = last(F.r[:, m]))
        aidx = argsample(F.a[:, m], N[m], scale = last(F.a[:, m]))
        aidx = Tuple.(CartesianIndices((Nf, Na))[aidx])
        idx = [(r = ri, f = fi, a = ai) for (ri, (fi, ai)) in zip(ridx, aidx)]
        push!(idxs, idx)
    end


    for (i, particle) in enumerate(cloud)
        m = ms[i]
        if !iszero(m)
            idx = pop!(idxs[m])
            value, state = birth(m, idx, ℓ, models, grid)
        else
            value, state = 1.0, State()
        end
        weights[i] *= F.Σm / value
        push!(particle, state)
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
