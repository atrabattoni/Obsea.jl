module Obsea

## Structures

export State, EmptyState, Trajectory, Metadata, Particle, Parameters
export move, logfk


struct State
    model::Int64
    frequency::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end


struct EmptyState end

Trajectory = Array{Union{State,EmptyState},1}


mutable struct Metadata
    weight::Float64
end


struct Particle
    trajectory::Trajectory
    metadata::Metadata
end


struct Parameters
    T::Float64
    q::Float64
    ps::Float64
    pb::Float64
end

## Motion Model

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

function move(state:EmptyState, params::Parameters)
    # TODO
end


function logfk(state::State, prevstate::State, params::Parameters)
    q, T, ps = params.q, params.T, params.ps
    dvx = state.vx - prevstate.vx
    dvy = state.vy - prevstate.vy
    log(ps) - (dvx^2 + dvy^2) / (q * T)^2  # TODO
end

function logfk(state::EmptyState, prevstate::State, params::Parameters)
    ps = params.ps
    log(1.0 - ps)
end

function logfk(state::State, prevstate::EmptyState, params::Parameters)
    pb = params.pb
    log(pb)
end

function logfk(state::EmptyState, prevstate::EmptyState, params::Parameters)
    pb = params.pb
    log(1.0 - pb)
end

end # module
