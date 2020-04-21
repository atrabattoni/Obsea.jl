@testset "types.jl" begin

    @testset "State" begin
        import Obsea: State, EmptyState, isempty
        state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
        @test state.model === 1
        @test state.frequency === 2.0
        @test state.x === 3.0
        @test state.y === 4.0
        @test state.vx === 5.0
        @test state.vy === 6.0
        ∅ = EmptyState()
        @test ∅.model === 0
        @test isempty(∅)
        @test !isempty(state)
    end

    @testset "Particle, Cloud and Weights" begin
        import Obsea: Particle, Cloud, Weights
        state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
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
        @test weights[1] == 1/10
        @test cloud[1] == Particle()
    end

    @testset "Parameters" begin
        import Obsea: Grid, Parameters
        grid = Grid(
            range(0.0, 1000.0, length = 5),
            range(0.0, 10.0, length = 7),
            range(0.0, 360.0, length = 9),
            (1:2),
        )
        params = Parameters(1.0, 0.1, 0.97, 0.03, 0.5, grid)
        @test params.T === 1.0
        @test params.q === 0.1
        @test params.ps === 0.97
        @test params.pb === 0.03
        @test params.pd === 0.5
    end

end
