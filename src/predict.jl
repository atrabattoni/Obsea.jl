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
        return move(state::State, params::Parameters)
    else
        return EmptyState()
    end
end

function transition(state::EmptyState, scan::Scan, params::Parameters)
    if rand() < params.pb
        return birth(scan::Scan, params::Parameters)
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
    return nothing
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
