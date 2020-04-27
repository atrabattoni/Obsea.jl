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

function final!(cloud)
    for particle in cloud
        popfirst!(particle)
    end
end

function estimate(cloud, Nm, Nt)
    m = zeros(Nm, Nt)
    x = zeros(Nt)
    y = zeros(Nt)
    n = zeros(Nt)
    for particle in cloud
        for (j, state) in enumerate(particle)
            if !isempty(state)
                n[j] += 1
                m[state.m, j] += 1
                x[j] += state.x
                y[j] += state.y
            end
        end
    end
    m ./= length(cloud)
    x ./= n
    y ./= n
    m, x, y
end
