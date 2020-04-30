function predict!(weights, cloud, prevcloud, ℓ, models, grid)
    mask = isempty.(cloud)
    @views birth!(weights[mask], cloud[mask], ℓ, models, grid)
    @views move!(cloud[.!mask], prevcloud[.!mask], models, grid)
    return cloud
end


## Birth

struct Index
    r::Int
    f::Int
    a::Int
    m::Int
end
Index() = Index(0, 0, 0, 0)

function birth!(weights, cloud, ℓ, models, grid)
    idxs = ℓ2idxs(ℓ, length(weights), grid)
    weights .*= idx2norm.(idxs, [ℓ])
    cloud .= idx2state.(idxs, [models], [grid])
end

function ℓ2idxs(ℓ, N, grid)
    @unpack Nr, Na, Nf, Nm = grid
    midxs = argsample(ℓ.m, N, scale = ℓ.Σm)
    Nc = counts(midxs, Nm)
    idxs = Vector{Index}(undef, N)
    j = 1
    for m = 1:Nm
        @views ridxs = argsample(grid.r .* ℓ.r[:, m], Nc[m])
        @views aidxs = argsample(vec(ℓ.a[:, :, m]), Nc[m])
        aidxs = Tuple.(CartesianIndices((Nf, Na))[aidxs])
        shuffle!(ridxs)
        shuffle!(aidxs)
        for (r, (f, a)) in zip(ridxs, aidxs)
            idxs[j] = Index(r, f, a, m)
            j += 1
        end
    end
    while j <= N
        idxs[j] = Index()
        j += 1
    end
    shuffle!(idxs)
end

function idx2norm(idx, ℓ)
    if !iszero(idx.m)
        return ℓ.Σm / ℓ.r[idx.r, idx.m] / ℓ.a[idx.f, idx.a, idx.m]
    else
        return ℓ.Σm
    end
end

function idx2state(idx, models, grid)
    if !iszero(idx.m)
        r = grid.r[idx.r]
        f = grid.f[idx.f]
        a = grid.a[idx.a]
        vr, va = randspeed(models[idx.m])
        x, y, = ra2xy(r, a)
        vx, vy = ra2xy(vr, va)
        return State(idx.m, f, x, y, vx, vy)
    else
        return State()
    end
end

function randspeed(model)
    @unpack vmin, vmax = model
    vr = vmin + (vmax - vmin) * rand()
    va = 2π * rand()
    va, vr
end


## Move

function move!(cloud, prevcloud, models, grid)
    Ps = [model.ps for model in models]
    Q = [model.q for model in models]
    T = grid.T
    kinematic!(cloud, prevcloud, T, Q)
    kill!(cloud, Ps)
    return cloud
end

function kinematic!(cloud, prevcloud, T, Q)
    dims = size(prevcloud)
    q = Q[prevcloud.m]
    ax = q .* randn(dims...)
    ay = q .* randn(dims...)
    cloud.x .= prevcloud.x .+ prevcloud.vx .* T .+ ax .* (T^2 / 2.0)
    cloud.y .= prevcloud.x .+ prevcloud.vy .* T .+ ay .* (T^2 / 2.0)
    cloud.vx .= prevcloud.vx .+ ax .* T
    cloud.vy .= prevcloud.vy .+ ay .* T
    return cloud
end

function kill!(prevcloud, Ps)
    dims = size(prevcloud)
    ps = Ps[prevcloud.m]
    mask = rand(dims...) .>= ps
    @views fill!(prevcloud[mask], State())
    return prevcloud
end


## TODO

# function logf(state, prevstate, models, grid)
#     @unpack T = grid
#     if isempty(state)
#         if isempty(prevstate)
#             return log(1.0 - pb)
#         else
#             return log(1.0 - ps)
#         end
#     else
#         if isempty(prevstate)
#             return log(pb)
#         else
#             @assert getmodel(state) == getmodel(prevstate)
#             @unpack q, pb, ps = models[getmodel(state)]
#             dvx = state.vx - prevstate.vx
#             dvy = state.vy - prevstate.vy
#             return log(ps) - (dvx^2 + dvy^2) / (q * T)^2
#         end
#     end
# end
