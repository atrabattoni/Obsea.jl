import Obsea: State, getmodel, isempty, Cloud, Weights, init, final!, estimate

@testset "particle" begin

    @testset "State" begin
        import Obsea: State, State, isempty
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
        @test isempty(∅)
        @test !isempty(state)
    end

    @testset "Particle, Cloud and Weights" begin
        import Obsea: Particle, Cloud, Weights
        state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = State()
        particle = [state, ∅]
        @test particle isa Particle
        cloud = [particle, particle]
        @test cloud isa Cloud
        weights = [0.1, 0.9]
        @test weights isa Weights
    end

    @testset "init" begin
        import Obsea.init
        weights, cloud = init(10)
        @test length(weights) == 10
        @test length(cloud) == 10
        @test weights[1] == 1 / 10
        @test cloud[1] == Particle([State()])
    end

    @testset "estimate" begin
        cloud = [
            [State(), State()],
            [State(1, 2.0, 2.0, 2.0, 2.0, 2.0), State(1, 3.0, 3.0, 3.0, 3.0, 3.0)],
            [State(1, 3.0, 3.0, 3.0, 3.0, 3.0), State(2, 3.0, 3.0, 3.0, 3.0, 3.0)],
        ]
        m = [
            (2 / 3) (1 / 3)
            0 (1 / 3)
        ]
        x = [2.5, 3.0]
        y = [2.5, 3.0]
        me, xe, ye = estimate(cloud, 2, 2)
        @test me ≈ m
        @test xe ≈ x
        @test ye ≈ y
    end

end
