@testset "particle.jl" begin

    @testset "ShipState" begin
        state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
        @test state.frequency === 2.0
        @test state.x === 3.0
        @test state.y === 4.0
        @test state.vx === 5.0
        @test state.vy === 6.0
    end

    @testset "Trajectory" begin
        state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        trajectory = Trajectory([state, ∅])
        @test trajectory == [state, ∅]
        @test trajectory[1] === state
    end

    @testset "Metadata" begin
        metadata = Metadata(0.1)
        @test metadata.weight === 0.1
    end

    @testset "Particle" begin
        state = ShipState(2.0, 3.0, 4.0, 5.0, 6.0)
        ∅ = EmptyState()
        trajectory = Trajectory([state, ∅])
        metadata = Metadata(0.1)
        particle = Particle(trajectory, metadata)
        @test particle.trajectory === trajectory
        @test particle.metadata === metadata
        particle.metadata.weight *= 2.0
        @test particle.metadata.weight === 0.2
    end

    @testset "Parameters" begin
        params = Parameters(1.0, 0.1, 0.97, 0.03, 0.5)
        @test params.T === 1.0
        @test params.q === 0.1
        @test params.ps === 0.97
        @test params.pb === 0.03
        @test params.pd === 0.5
    end

end
