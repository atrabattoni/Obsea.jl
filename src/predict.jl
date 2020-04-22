function predict!(cloud, scan, model)
    for particle in cloud
        push!(particle, transition(particle[end], scan, model))
    end
end


function transition(state, scan, model)
    @unpack pb, ps = model
    if isempty(state)
        if rand() < pb
            return birth(scan, model)
        else
            return EmptyState()
        end
    else
        if rand() < ps
            return move(state, model)
        else
            return EmptyState()
        end
    end
end


function move(state, model)
    @unpack q, T = model
    ax = q * randn()
    ay = q * randn()
    model = state.model
    frequency = state.frequency
    x = state.x + state.vx * T + ax * T^2 / 2.0
    y = state.y + state.vy * T + ay * T^2 / 2.0
    vx = state.vx + ax * T
    vy = state.vy + ay * T
    State(model, frequency, x, y, vx, vy)
end


function birth(scan, model)
    @unpack grid = model
    idx_r = searchsortedfirst(scan.cdf_r, rand(), lt = <=)
    idx_fam = searchsortedfirst(scan.cdf_fam, rand(), lt = <=)
    cidx_fam = CartesianIndex((
        length(grid.range_f),
        length(grid.range_a),
        length(grid.range_m),
    ))
    r = grid.range_r[idx_r]
    f = grid.range_f[cidx_fam[1]]
    a = grid.range_a[cidx_fam[2]]
    m = cidx_fam[3]
    x = r * sin(a)
    y = r * cos(a)
    vx = 10.0 * rand()
    vy = 10.0 * rand()
    State(m, f, x, y, vx, vy)
end


function logf(state, prevstate, model)
    @unpack pb, ps, q, T = model
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
            dvx = state.vx - prevstate.vx
            dvy = state.vy - prevstate.vy
            return log(ps) - (dvx^2 + dvy^2) / (q * T)^2  # TODO
        end
    end
end
