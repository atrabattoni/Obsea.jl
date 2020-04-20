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

const AnyState = Union{EmptyState,ShipState,WhaleState}

const Trajectory = Vector{AnyState}

mutable struct Metadata
    weight::Float64
end

struct Particle
    trajectory::Trajectory
    metadata::Metadata
end

const Cloud = Vector{Particle}

struct Grid
    range_r::AbstractRange
    range_f::AbstractRange
    range_a::AbstractRange
    range_m::UnitRange{Int64}
end

struct Scan
    cdf_r::Array{Float64,1}
    cdf_fam::Array{Float64,1}
    itp_r::ScaledInterpolation
    itp_fam::ScaledInterpolation
end

@with_kw struct Parameters
    T::Float64 = NaN
    q::Float64 = NaN
    ps::Float64 = NaN
    pb::Float64 = NaN
    pd::Float64 = NaN
    grid::Grid = Grid(0.0:0.0, 0.0:0.0, 0.0:0.0, 0:0)
end

function init(N)
    cloud = Cloud([Particle(Trajectory(), Metadata(1 / N)) for i = 1:N])
end
