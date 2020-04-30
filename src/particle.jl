struct State
    m::Int64
    f::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end
State() = State(0, 0.0, 0.0, 0.0, 0.0, 0.0)
getmodel(state::State) = state.m
isempty(state::State) = iszero(getmodel(state))

function init(Nt, Np)
    weights = fill(1 / Np, Np)
    particles = StructArray{State}(undef, Nt, Np)
    fill!(particles, State())
    return weights, particles
end

function estimate(particles)
    n = sum(.!iszero.(particles.m))
    x = sum(particles.x, dims = 2) / n
    y = sum(particles.y, dims = 2) / n
    return vec(x), vec(y)
end
