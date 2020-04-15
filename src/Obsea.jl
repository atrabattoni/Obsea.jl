module Obsea

export State, Trajectory, Metadata, Particle

struct State
    model::Int64
    frequency::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end

Trajectory = Array{State, 1}

mutable struct Metadata
    weight::Float64
end

struct Particle
    trajectory::Trajectory
    metadata::Metadata
end

end # module
