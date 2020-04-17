@testset "predict.jl" begin

    state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
    movedstate = ShipState(2.0, 8.0, 10.0, 5.0, 6.0)
    ∅ = EmptyState()
    grid = Grid(
        range(0.0, 1000.0, length = 5),
        range(0.0, 10.0, length = 7),
        range(0.0, 360.0, length = 9),
        (1:2),
    )
    scan = Scan(ones(5), ones(7, 9, 2), grid)

    @testset "move" begin
        params = Parameters(1.0, 0.0, 0.97, 0.03, 0.5, grid)
        @test move(state, params) == movedstate
    end

    @testset "transition" begin
        params = Parameters(1.0, 0.0, 0.0, 0.0, 1.0, grid)
        @test transition(state, scan, params) == EmptyState()
        @test transition(∅, scan, params) == EmptyState()
        params = Parameters(1.0, 0.0, 1.0, 1.0, 1.0, grid)
        @test transition(state, scan, params) == movedstate
        @test typeof(transition(∅, scan, params)) <: State
    end

    @testset "logf" begin
        params = Parameters(1.0, 0.1, 0.97, 0.03, 0.5, grid)
        @test logf(state, state, params) === log(params.ps)  # TODO: diff state
        @test logf(∅, state, params) === log(1.0 - params.ps)
        @test logf(state, ∅, params) === log(params.pb)
        @test logf(∅, ∅, params) === log(1.0 - params.pb)
    end

    @testset "predict" begin
        metadata = Metadata(0.1)
        trajectory = Trajectory([state])
        particle = Particle(trajectory, metadata)
        cloud = Cloud([particle])

        params = Parameters(1.0, 0.0, 1.0, 0.0, 0.5, grid) # TODO: pb ≠ 0
        predict!(cloud, scan, params)
        @test length(particle.trajectory) === 2
        @test particle.trajectory[2] == movedstate

        params = Parameters(1.0, 0.0, 0.0, 0.0, 0.5, grid)
        predict!(cloud, scan, params)
        @test length(particle.trajectory) === 3
        @test particle.trajectory[3] == ∅
    end

end
