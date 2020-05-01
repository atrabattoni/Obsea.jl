import Obsea:
    predict!,
    Index,
    birth!,
    ℓ2idxs,
    idx2norm,
    idx2state,
    randspeed,
    move!,
    kinematic!,
    survive

@testset "predict.jl" begin
    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    @unpack T, Nr, Nf, Na, Nm = grid
    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    movedstate = State(1, 2.0, 3.0 + 5.0 * T, 4.0 + 6.0 * T, 5.0, 6.0)
    ∅ = State()
    life = Model(
        name = "life",
        q = 0.0,
        vmin = 0.0,
        vmax = 0.0,
        ps = 1.0,
        pb = 1.0,
        pd = NaN,
        lam = NaN,
        mrl = NaN,
        n = 1,
    )
    death = Model(
        name = "death",
        q = 0.0,
        vmin = 0.0,
        vmax = 0.0,
        ps = 0.0,
        pb = 0.0,
        pd = NaN,
        lam = NaN,
        mrl = NaN,
        n = 1,
    )
    lifedeath = StructVector([life, death])
    deathlife = StructVector([death, life])
    deathdeath = StructVector([death, death])

    @testset "move" begin
        prevcloud = StructVector([state])
        cloud = StructVector([State()])
        move!(cloud, prevcloud, lifedeath, grid)
        @test cloud == StructVector([movedstate])
        cloud = StructVector([State()])
        move!(cloud, prevcloud, deathdeath, grid)
        @test cloud.m == [0]
    end

    @testset "birth!, ℓ2idxs, idx2state, idx2norm" begin
        idx = Index(57, 33, 101, 2)
        value = 10.0
        N = 3
        ℓ = (
            r = zeros(Nr, Nm),
            a = zeros(Nf, Na, Nm),
            m = zeros(Nm),
            Σm = value * value / Nr / Nf / Na,
        )
        ℓ.r[idx.r, idx.m] = value
        ℓ.a[idx.f, idx.a, idx.m] = value
        ℓ.m[idx.m] = value * value / Nr / Nf / Na
        r = grid.r[idx.r]
        f = grid.f[idx.f]
        a = grid.a[idx.a]
        m = idx.m
        x = r * sin(a)
        y = r * cos(a)
        # ℓ2idxs
        idxs = ℓ2idxs(ℓ, N, grid)
        @test length(idxs) == N
        @test idxs == fill(idx, N)
        # idx2state
        s = idx2state(idxs[1], deathlife, grid)
        @test s.m == m
        @test s.f == f
        @test s.x == x
        @test s.y == y
        @test s.vx == 0.0
        @test s.vy == 0.0

        # idx2norm
        norm = idx2norm(idxs[1], ℓ)
        @test norm == ℓ.Σm / ℓ.r[idx.r, idx.m] / ℓ.a[idx.f, idx.a, idx.m]
        # birth
        weights = [1.0]
        cloud = StructVector([State()])
        birth!(weights, cloud, ℓ, deathlife, grid)
        @test weights == [ℓ.Σm / ℓ.r[idx.r, idx.m] / ℓ.a[idx.f, idx.a, idx.m]]
        s = cloud[1]
        @test s.m == m
        @test s.f == f
        @test s.x == x
        @test s.y == y
        @test s.vx == 0.0
        @test s.vy == 0.0
    end

    @testset "isdead, survive, predict" begin
        ℓ = Likelihood(ones(Nr, Nm), ones(Nf, Na, Nm), ones(Nm), Nm)
        prevcloud = StructVector([state])
        cloud = StructVector([State()])
        weights = [1.0]

        @test isdead.(prevcloud) == [false]
        @test survive(prevcloud, [1.0]) == [true]
        move!(cloud, prevcloud, lifedeath, grid)
        @test cloud[1] == movedstate
        predict!(weights, cloud, prevcloud, ℓ, lifedeath, grid)
        @test cloud[1] == movedstate

        prevcloud = cloud
        cloud = StructVector([State()])
        @test isdead.(prevcloud) == [false]
        @test survive(prevcloud, [0.0]) == [false]
        predict!(weights, cloud, prevcloud, ℓ, deathdeath, grid)
        @test isdead.(cloud) == [true]
    end

end
