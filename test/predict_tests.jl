@testset "predict.jl" begin

    import Obsea.Scan
    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    movedstate = State(1, 2.0, 8.0, 10.0, 5.0, 6.0)
    ∅ = EmptyState()
    grid = Grid(
        range(0.0, 1000.0, length = 5),
        range(0.0, 10.0, length = 7),
        range(0.0, 360.0, length = 9),
        (1:2),
    )
    scan = Scan(ones(5), ones(7, 9, 2), grid)

    @testset "move" begin
        import Obsea.move
        model = Model(1.0, 0.0, 0.97, 0.03, 0.5, grid)
        @test move(state, model) == movedstate
    end

    @testset "transition" begin
        import Obsea.transition
        model = Model(1.0, 0.0, 0.0, 0.0, 1.0, grid)
        @test isempty(transition(state, scan, model))
        @test isempty(transition(∅, scan, model))
        model = Model(1.0, 0.0, 1.0, 1.0, 1.0, grid)
        @test transition(state, scan, model) == movedstate
        @test transition(∅, scan, model) isa State
    end

    @testset "logf" begin
        import Obsea.logf
        model = Model(1.0, 0.1, 0.97, 0.03, 0.5, grid)
        @test logf(state, state, model) === log(model.ps)  # TODO: diff state
        @test logf(∅, state, model) === log(1.0 - model.ps)
        @test logf(state, ∅, model) === log(model.pb)
        @test logf(∅, ∅, model) === log(1.0 - model.pb)
    end

    @testset "predict" begin
        import Obsea.predict!
        particle = Particle([state])
        cloud = Cloud([particle])

        model = Model(1.0, 0.0, 1.0, 0.0, 0.5, grid) # TODO: pb ≠ 0
        predict!(cloud, scan, model)
        @test length(particle) === 2
        @test particle[2] == movedstate

        model = Model(1.0, 0.0, 0.0, 0.0, 0.5, grid)
        predict!(cloud, scan, model)
        @test length(particle) === 3
        @test isempty(particle[3])
    end

end
