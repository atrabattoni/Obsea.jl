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


struct Grid
    rrange::AbstractRange
    arange::AbstractRange
    frange::AbstractRange  
end


struct Scan
    yr::Array{Float64,1}
    ya::Array{Float64,3}
end
