import Obsea: State, getmodel, isdead, init, estimate

@testset "particle" begin

    @testset "State, isdead, getmodel" begin
        state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
        @test state.m === 1
        @test state.f === 2.0
        @test state.x === 3.0
        @test state.y === 4.0
        @test state.vx === 5.0
        @test state.vy === 6.0
        ∅ = State()
        @test getmodel(∅) === 0
        @test getmodel(state) === 1
        @test isdead(∅)
        @test !isdead(state)
        @test isdead.(StructVector([state, State()])) == [false, true]
    end

    @testset "init" begin
        weights, particles = init(10, 3)
        @test length(weights) == 3
        @test size(particles) == (10, 3)
        @test weights[1] == 1 / 3
        @test particles[1, 1] == State()
    end

    @testset "estimate" begin
        sa = State(1, 2.0, 2.0, 2.0, 2.0, 2.0)
        sb = State(1, 3.0, 3.0, 3.0, 3.0, 3.0)
        sc = State(2, 3.0, 3.0, 3.0, 3.0, 3.0)
        particles = StructArray([
            State() sa sb
            State() sb sc
        ])
        x = [2.5, 3.0]
        y = [2.5, 3.0]
        xe, ye = estimate(particles)
        @test xe ≈ x
        @test ye ≈ y
    end

end
