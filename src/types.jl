struct State
    model::Int64
    frequency::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end

EmptyState() = State(0, NaN, NaN, NaN, NaN, NaN)
isempty(s::State) = (s.model == 0)

const Particle = Vector{State}
const Cloud = Vector{Particle}
const Weights = Vector{Float64}

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

struct Model
    T::Float64
    q::Float64
    ps::Float64
    pb::Float64
    pd::Float64
    grid::Grid
end

function init(N)
    weights = Weights(fill(1 / N, N))
    cloud = Cloud([Particle() for i = 1:N])
    return (weights, cloud)
end
