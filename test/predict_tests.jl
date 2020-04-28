import Obsea: predict!, transition, move, birth, logf

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

    @testset "birth" begin
        ℓ = (r = ones(Nr, Nm), a = ones(Nf, Na, Nm), m = ones(Nm))

        F = distribution(ℓ, [life, death], grid)
        @test getmodel(last(birth(ℓ, F, [life, death], grid))) === 1

        F = distribution(ℓ, [death, life], grid)
        @test getmodel(last(birth(ℓ, F, [death, life], grid))) === 2

        F = distribution(ℓ, [death, death], grid)
        @test isempty(last(birth(ℓ, F, [death, death], grid)))

        F = distribution(ℓ, [half, half], grid)
        @test !isempty(last(birth(ℓ, F, [half, half], grid)))

        F = distribution(ℓ, [life, death], grid)
        @test first(birth(ℓ, F, [life, death], grid)) ≈ 1.0

        F = distribution(ℓ, [death, life], grid)
        @test first(birth(ℓ, F, [death, life], grid)) ≈ 1.0

        F = distribution(ℓ, [death, death], grid)
        @test first(birth(ℓ, F, [death, death], grid)) ≈ 1.0

        F = distribution(ℓ, [half, half], grid)
        @test first(birth(ℓ, F, [half, half], grid)) ≈ 1.0

        F = distribution(ℓ, [half, death], grid)
        @test first(birth(ℓ, F, [half, death], grid)) ≈ 1.0

        F = distribution(ℓ, [death, half], grid)
        @test first(birth(ℓ, F, [death, half], grid)) ≈ 1.0

        F = distribution(ℓ, [death, death], grid)
        @test first(birth(ℓ, F, [death, death], grid)) ≈ 1.0


        ℓ = (r = 2 * ones(Nr, Nm), a = 2 * ones(Nf, Na, Nm), m = 2 * ones(Nm))
        normalization = 1
        while true
            F = distribution(ℓ, [half, death], grid)
            normalization, s = birth(ℓ, F, [half, death], grid)
            if isempty(s)
                break
            end
        end
        @test normalization > 1.0
        while true
            F = distribution(ℓ, [half, death], grid)
            normalization, s = birth(ℓ, F, [half, death], grid)
            if !isempty(s)
                break
            end
        end
        @test normalization < 1.0

        ℓ = (r = ones(Nr, Nm) / 2, a = ones(Nf, Na, Nm) / 2, m = ones(Nm) / 2)
        while true
            F = distribution(ℓ, [half, death], grid)
            normalization, s = birth(ℓ, F, [half, death], grid)
            if isempty(s)
                break
            end
        end
        @test normalization < 1.0
        while true
            F = distribution(ℓ, [half, death], grid)
            normalization, s = birth(ℓ, F, [half, death], grid)
            if !isempty(s)
                break
            end
        end
        @test normalization > 1.0

        ℓ = (r = zeros(Nr, Nm), a = zeros(Nf, Na, Nm), m = zeros(Nm))
        ℓ.r[101, 2] = 10.0
        ℓ.a[101, 101, 2] = 10.0
        ℓ.m[2] = 10.0
        r = grid.r[101]
        a = grid.a[101]
        F = distribution(ℓ, [death, life], grid)
        _, b = birth(ℓ, F, [death, life], grid)
        @test b.f == grid.f[101]
        @test b.x == r * sin(a)
        @test b.y == r * cos(a)

    end

    # @testset "transition" begin
    #     ℓ = (r = ones(Nr, Nm), a = ones(Nf, Na, Nm), m = ones(Nm))
    #
    #     F = distribution(ℓ, [life, death], grid)
    #     @test last(transition(state, ℓ, F, [life, death], grid)) == movedstate
    #
    #     F = distribution(ℓ, [life, death], grid)
    #     @test getmodel(last(transition(∅, ℓ, F, [life, death], grid))) == 1
    #
    #     F = distribution(ℓ, [death, life], grid)
    #     @test getmodel(last(transition(∅, ℓ, F, [death, life], grid))) == 2
    #
    #     F = distribution(ℓ, [death, death], grid)
    #     @test isempty(last(transition(state, ℓ, F, [death, death], grid)))
    #
    #     F = distribution(ℓ, [death, death], grid)
    #     @test isempty(last(transition(∅, ℓ, F, [death, death], grid)))
    # end

    @testset "predict" begin
        ℓ = (r = ones(Nr, Nm), a = ones(Nf, Na, Nm), m = ones(Nm))
        particle = [state]
        cloud = [particle]
        weights = [1.0]
        F = distribution(ℓ, [life, death], grid)
        predict!(weights, cloud, ℓ, F, [life, death], grid)
        @test length(particle) === 2
        @test particle[2] == movedstate

        F = distribution(ℓ, [death, death], grid)
        predict!(weights, cloud, ℓ, F, [death, death], grid)
        @test length(particle) === 3
        @test isempty(particle[3])
    end
    #
    # @testset "logf" begin
    #     model = Model(1.0, 0.1, 0.97, 0.03, 0.5, grid)
    #     @test logf(state, state, model) === log(model.ps)  # TODO: diff state
    #     @test logf(∅, state, model) === log(1.0 - model.ps)
    #     @test logf(state, ∅, model) === log(model.pb)
    #     @test logf(∅, ∅, model) === log(1.0 - model.pb)
    # end


end
