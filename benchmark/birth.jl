using StructArrays
using BenchmarkTools
using Profile

struct Index
    r::Int
    f::Int
    a::Int
    m::Int
end
Index() = Index(0, 0, 0, 0)

struct Likelihood
    r::Array{Float64,2}
    a::Array{Float64,3}
    m::Vector{Float64}
    Σm::Float64
end

function idx2norm(idx, ℓ)
    if !iszero(idx.m)
        return ℓ.Σm / ℓ.r[idx.r, idx.m] / ℓ.a[idx.f, idx.a, idx.m]
    else
        return ℓ.Σm
    end
end



Nr = 513
Nf = 266
Na = 361
Nm = 2
Np = 1000

ℓ = Likelihood(ones(Nr, Nm), ones(Nf, Na, Nm), ones(Nm), 1.0 + Nm)

idx = StructVector{Index}((
    rand(1:Nr, Np),
    rand(1:Nf, Np),
    rand(1:Na, Np),
    rand(1:Nm, Np),
))

@btime idx2norm.(idx, [ℓ])
@btime for i = 1:Np idx2norm(idx[i], ℓ) end



function birth!(mask, weights, cloud, ℓ, models, grid)
    idxs = ℓ2idxs(ℓ, count(mask), grid)
    weights[mask] .*= idx2norm.(idxs, ℓ)
    cloud[mask] .= idx2state.(idxs, models, grid)
end

function ℓ2idxs(ℓ, Np, grid)
    @unpack Nr, Na, Nf, Nm = grid
    idxs = StructVector{Index}(undef, Np)
    idxs.m .= argsample(ℓ.m, Np, scale = ℓ.Σm)
    Nc = counts(idxs.m, Nm)
    Struct
    j = 1
    for m = 1:Nm
        @views idxs.r[j:j+N[m]-1] =
            shuffle!(argsample(grid.r .* ℓ.r[:, m], Nc[m]))
        @views aidxs = argsample(vec(ℓ.a[:, :, m]), Nc[m])
        aidxs = Tuple.(CartesianIndices((Nf, Na))[aidxs])
        shuffle!(ridxs)
        shuffle!(aidxs)

    end
    while j <= Np
        idxs[j] = Index()
        j += 1
    end
    shuffle!(idxs)
end

function counts(x, Np)
    d = Dict(i => 0 for i = 0:Np)
    for xi in x
        d[xi] += 1
    end
    d
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
