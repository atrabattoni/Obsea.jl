struct State
    m::Int64
    f::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end

State() = State(0, NaN, NaN, NaN, NaN, NaN)
getmodel(s) = s.m
Base.isempty(s::State) = iszero(getmodel(s))

const Particle = Vector{State}
const Cloud = Vector{Particle}
const Weights = Vector{Float64}

function init(N)
    weights = Weights(fill(1 / N, N))
    cloud = Cloud([Particle([State()]) for i = 1:N])
    return (weights, cloud)
end
