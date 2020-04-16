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
    model = state.model
    frequency = state.frequency
    x = state.x + state.vx * T + ax * T^2 / 2.0
    y = state.y + state.vy * T + ay * T^2 / 2.0
    vx = state.vx + ax * T
    vy = state.vy + ay * T
    State(model, frequency, x, y, vx, vy)
end


function birth(scan::Scan, params::Parameters)
    return nothing
end
