@testset "predict.jl" begin

    @testset "move" begin
        state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        params = Parameters(1.0, 0.0, 0.97, 0.03, 0.5)
        @test move(state, params) == State(1, 2.0, 8.0, 10.0, 5.0, 6.0)
    end

    @testset "logf" begin
        state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        params = Parameters(1.0, 0.1, 0.97, 0.03, 0.5)
        state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        @test logf(state, state, params) === log(params.ps)  # TODO
        @test logf(∅, state, params) === log(1.0 - params.ps)
        @test logf(state, ∅, params) === log(params.pb)
        @test logf(∅, ∅, params) === log(1.0 - params.pb)
    end

end
