function predict!(cloud::Cloud, scan::Scan, params::Parameters)
    for particle ∈ cloud
        push!(
            particle.trajectory,
            transition(particle.trajectory[end], scan, params),
        )
    end
end


function transition(state::State, scan::Scan, params::Parameters)
    if rand() < params.ps
        return move(state, params)
    else
        return EmptyState()
    end
end

function transition(state::EmptyState, scan::Scan, params::Parameters)
    if rand() < params.pb
        return birth(scan, params)
    else
        return EmptyState()
    end
end


function move(state::State, params::Parameters)
    q, T = params.q, params.T
    ax = q * randn()
    ay = q * randn()
    frequency = state.frequency
    x = state.x + state.vx * T + ax * T^2 / 2.0
    y = state.y + state.vy * T + ay * T^2 / 2.0
    vx = state.vx + ax * T
    vy = state.vy + ay * T
    typeof(state)(frequency, x, y, vx, vy)
end


function birth(scan::Scan, params::Parameters)
    grid = params.grid
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
    if m == 1
        return ShipState(f, x, y, vx, vy)
    elseif m == 2
        return WhaleState(f, x, y, vx, vy)
    end
end


function logf(state::State, prevstate::State, params::Parameters)
    q, T, ps = params.q, params.T, params.ps
    dvx = state.vx - prevstate.vx
    dvy = state.vy - prevstate.vy
    log(ps) - (dvx^2 + dvy^2) / (q * T)^2  # TODO
end

function logf(state::EmptyState, prevstate::State, params::Parameters)
    ps = params.ps
    log(1.0 - ps)
end

function logf(state::State, prevstate::EmptyState, params::Parameters)
    pb = params.pb
    log(pb)
end

function logf(state::EmptyState, prevstate::EmptyState, params::Parameters)
    pb = params.pb
    log(1.0 - pb)
end
