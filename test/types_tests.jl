@testset "particle.jl" begin

    state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
    ∅ = EmptyState()
    trajectory = Trajectory([state, ∅])
    metadata = Metadata(0.1)
    particle = Particle(trajectory, metadata)
    cloud = Cloud([particle, particle])
    params = Parameters(1.0, 0.1, 0.97, 0.03, 0.5)

    @testset "ShipState" begin
        @test state.frequency === 2.0
        @test state.x === 3.0
        @test state.y === 4.0
        @test state.vx === 5.0
        @test state.vy === 6.0
    end

    @testset "Trajectory" begin
        @test trajectory == [state, ∅]
        @test trajectory[1] === state
    end

    @testset "Metadata" begin
        @test metadata.weight === 0.1
    end

    @testset "Particle" begin
        @test particle.trajectory === trajectory
        @test particle.metadata === metadata
        particle.metadata.weight *= 2.0
        @test particle.metadata.weight === 0.2
    end

    @testset "Cloud" begin
        @test cloud[1] == particle
    end

    @testset "init" begin
        cloud = init(10)
        @test length(cloud) === 10
        @test cloud[1].trajectory == Trajectory()
        @test cloud[1].metadata.weight == 1/10
    end

    @testset "Parameters" begin
        @test params.T === 1.0
        @test params.q === 0.1
        @test params.ps === 0.97
        @test params.pb === 0.03
        @test params.pd === 0.5
    end

end
