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
