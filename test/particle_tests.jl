import Obsea: State, getmodel, isempty, Cloud, Weights, init

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
        @test cloud[1] == Particle()
    end

end
