module Obsea

## Structures

export State, EmptyState, Trajectory, Metadata, Particle, Parameters
export qk, logfk

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

function qk(state::State, params::Parameters)
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

function logfk(state::State, prevstate::State, params::Parameters)
    q, T = params.q, params.T
    dvx = state.vx - prevstate.vx
    dvy = state.vy - prevstate.vy
    -(dvx^2 + dvy^2) / (q * T)^2  # TODO
end

end # module
