module Obsea

export State, Particle

struct State
    model::Int64
    frequency::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end

struct Particle
    trajectory::Array{State, 1}
    logLikelihood::Float64
    weight::Float64
end

end # module
