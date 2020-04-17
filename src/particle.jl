abstract type State end

struct ShipState <: State
    frequency::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end

struct WhaleState <: State
    frequency::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end

struct EmptyState end

AnyState = Union{EmptyState,ShipState,WhaleState}

Trajectory = Vector{AnyState}

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
    pd::Float64
end
