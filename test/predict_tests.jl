import Obsea:
    predict!,
    birth!,
    ℓ2idxs,
    counts,
    idx2state,
    randspeed,
    idx2norm,
    move!,
    move,
    logf

@testset "predict.jl" begin

    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    @unpack T = grid
    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    movedstate = State(1, 2.0, 3.0 + 5.0 * T, 4.0 + 6.0 * T, 5.0, 6.0)
    ∅ = State()
    life = Model(
        name = "life",
        q = 0.0,
        vmin = 0.0,
        vmax = 1.0,
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
        vmax = 1.0,
        ps = 0.0,
        pb = 0.0,
        pd = NaN,
        lam = NaN,
        mrl = NaN,
        n = 1,
    )
    half = Model(
        name = "death",
        q = 0.0,
        vmin = 0.0,
        vmax = 1.0,
        ps = 0.5,
        pb = 0.5,
        pd = NaN,
        lam = NaN,
        mrl = NaN,
        n = 1,
    )
    Nr, Nf, Na, Nm = (grid.Nr, grid.Nf, grid.Na, length(1:grid.Nm))

    @testset "move" begin
        @test move(state, [life, death], grid) == movedstate
        @test isempty(move(state, [death, death], grid))
    end

    @testset "distribution, ℓ2idxs, idx2state, idx2norm" begin
        idx = (r = 57, f = 33, a = 101, m = 2)
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
        s = idx2state(idxs[1], [death, life], grid)
        @test s.m == m
        @test s.f == f
        @test s.x == x
        @test s.y == y
        # idx2norm
        norm = idx2norm(idxs[1], ℓ)
        @test norm == ℓ.Σm / ℓ.r[idx.r, idx.m] / ℓ.a[idx.f, idx.a, idx.m]
    end

    @testset "predict" begin
        ℓ = (r = ones(Nr, Nm), a = ones(Nf, Na, Nm), m = ones(Nm), Σm = Nm)
        particle = [state]
        cloud = [particle]
        weights = [1.0]
        predict!(weights, cloud, ℓ, [life, death], grid)
        @test length(particle) === 2
        @test particle[2] == movedstate

        predict!(weights, cloud, ℓ, [death, death], grid)
        @test length(particle) === 3
        @test isempty(particle[3])
    end

end
