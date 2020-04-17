@testset "predict.jl" begin

    @testset "move" begin
        state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        params = Parameters(1.0, 0.0, 0.97, 0.03, 0.5)
        @test move(state, params) == ShipState(2.0, 8.0, 10.0, 5.0, 6.0)
    end

    @testset "transition" begin
        state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        grid = Grid(
            range(0.0, 1000.0, length = 5),
            range(0.0, 10.0, length = 7),
            range(0.0, 360.0, length = 9),
            (1:2),
        )
        scan = Scan(ones(5), ones(7, 9, 2), grid)
        params = Parameters(1.0, 0.0, 0.0, 0.0, 1.0)
        @test transition(state, scan, params) == EmptyState()
        @test transition(∅, scan, params) == EmptyState()
        params = Parameters(1.0, 0.0, 1.0, 1.0, 1.0)
        @test transition(state, scan, params) ==
              ShipState(2.0, 8.0, 10.0, 5.0, 6.0)
        @test transition(∅, scan, params) == nothing
    end

    @testset "logf" begin
        state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        params = Parameters(1.0, 0.1, 0.97, 0.03, 0.5)
        state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        @test logf(state, state, params) === log(params.ps)  # TODO
        @test logf(∅, state, params) === log(1.0 - params.ps)
        @test logf(state, ∅, params) === log(params.pb)
        @test logf(∅, ∅, params) === log(1.0 - params.pb)
    end

end
