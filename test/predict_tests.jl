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
        vmax = Inf,
        ps = 1.0,
        pb = 1.0,
        pd = NaN,
        lam = NaN,
        mrl = NaN,
        n = NaN,
    )
    death = Model(
        name = "death",
        q = 0.0,
        vmin = 0.0,
        vmax = Inf,
        ps = 0.0,
        pb = 0.0,
        pd = NaN,
        lam = NaN,
        mrl = NaN,
        n = NaN,
    )
    half = Model(
        name = "death",
        q = 0.0,
        vmin = 0.0,
        vmax = Inf,
        ps = 0.5,
        pb = 0.5,
        pd = NaN,
        lam = NaN,
        mrl = NaN,
        n = NaN,
    )
    dims = (
        length(grid.rrange),
        length(grid.frange),
        length(grid.arange),
        length(grid.mrange),
    )

    @testset "move" begin
        @test last(move(state, [life, death], grid)) == movedstate
        @test isempty(last(move(state, [death, death], grid)))
    end

    @testset "birth" begin
        ℓ = ones(dims)
        @test getmodel(last(birth(ℓ, [life, death], grid))) === 1
        @test getmodel(last(birth(ℓ, [death, life], grid))) === 2
        @test isempty(last(birth(ℓ, [death, death], grid)))
        @test !isempty(last(birth(ℓ, [half, half], grid)))
        @test first(birth(ℓ, [life, death], grid)) ≈ 1.0
        @test first(birth(ℓ, [death, life], grid)) ≈ 1.0
        @test first(birth(ℓ, [death, death], grid)) ≈ 1.0
        @test first(birth(ℓ, [half, half], grid)) ≈ 1.0
        @test first(birth(ℓ, [half, death], grid)) ≈ 1.0
        @test first(birth(ℓ, [death, half], grid)) ≈ 1.0
        @test first(birth(ℓ, [death, death], grid)) ≈ 1.0


        ℓ = 2 * ones(dims)
        normalization = 1
        while true
            normalization, s = birth(ℓ, [half, death], grid)
            if isempty(s)
                break
            end
        end
        @test normalization > 1.0
        while true
            normalization, s = birth(ℓ, [half, death], grid)
            if !isempty(s)
                break
            end
        end
        @test normalization < 1.0

        ℓ = ones(dims) / 2
        while true
            normalization, s = birth(ℓ, [half, death], grid)
            if isempty(s)
                break
            end
        end
        @test normalization < 1.0
        while true
            normalization, s = birth(ℓ, [half, death], grid)
            if !isempty(s)
                break
            end
        end
        @test normalization > 1.0

    end

    @testset "transition" begin
        ℓ = ones(dims)
        @test last(transition(state, ℓ, [life, death], grid)) == movedstate
        @test getmodel(last(transition(∅, ℓ, [life, death], grid))) == 1
        @test getmodel(last(transition(∅, ℓ, [death, life], grid))) == 2
        @test isempty(last(transition(state, ℓ, [death, death], grid)))
        @test isempty(last(transition(∅, ℓ, [death, death], grid)))
    end

    @testset "predict" begin
        ℓ = ones(dims)
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
    #
    # @testset "logf" begin
    #     model = Model(1.0, 0.1, 0.97, 0.03, 0.5, grid)
    #     @test logf(state, state, model) === log(model.ps)  # TODO: diff state
    #     @test logf(∅, state, model) === log(1.0 - model.ps)
    #     @test logf(state, ∅, model) === log(model.pb)
    #     @test logf(∅, ∅, model) === log(1.0 - model.pb)
    # end


end
