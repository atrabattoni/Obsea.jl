function predict!(weights, cloud, ℓ, models, grid)
    mask = isempty(cloud)
    birth!(view(weights, mask), view(cloud, mask), ℓ, models, grid)
    move!(view(cloud, .!mask), models, grid)
end

function birth!(weights, cloud, ℓ, models, grid)
    idxs = ℓ2idxs(ℓ, length(cloud), grid)
    for (j, idx, particle) in zip(1:length(weights), idxs, cloud)
        weights[j] *= idx2norm(idx, ℓ)
        push!(particle, idx2state(idx, models, grid))
    end
end

function ℓ2idxs(ℓ, N, grid)
    @unpack Nr, Na, Nf, Nm = grid
    midxs = argsample(ℓ.m, N, scale = ℓ.Σm)
    Nc = counts(midxs, Nm)
    idxs = Array{NamedTuple{(:r, :f, :a, :m),NTuple{4,Int}}}(undef, N)
    j = 1
    for m = 1:Nm
        @views ridxs = argsample(ℓ.r[:, m], Nc[m])
        @views faidxs = argsample(vec(ℓ.a[:, :, m]), Nc[m])
        faidxs = Tuple.(CartesianIndices((Nf, Na))[faidxs])
        for (r, (f, a)) in zip(ridxs, faidxs)
            idxs[j] = (r = r, f = f, a = a, m = m)
            j += 1
        end
    end
    while j <= N
        idxs[j] = (r = 0, f = 0, a = 0, m = 0)
        j += 1
    end
    shuffle!(idxs)
end

function counts(x, N)
    d = Dict(i => 0 for i = 0:N)
    for xi in x
        d[xi] += 1
    end
    d
end

function idx2state(idx, models, grid)
    if !iszero(idx.m)
        r = grid.r[idx.r]
        f = grid.f[idx.f]
        a = grid.a[idx.a]
        vr, va = randspeed(models[idx.m])
        x, y, = ra2xy(r, a)
        vx, vy = ra2xy(vr, va)
        return State(idx.m, f, x, y, vx, vy)
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

function idx2norm(idx, ℓ)
    if !iszero(idx.m)
        return ℓ.Σm / ℓ.r[idx.r, idx.m] / ℓ.a[idx.f, idx.a, idx.m]
    else
        return ℓ.Σm
    end
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
