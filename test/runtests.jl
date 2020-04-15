using Obsea
using Test

@testset "Obsea.jl" begin

@testset "Structures" begin
    params =
    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    @testset "State" begin
        @test state.model === 1
        @test state.frequency === 2.0
        @test state.x === 3.0
        @test state.y === 4.0
        @test state.vx === 5.0
        @test state.vy === 6.0
    end
    ∅ = EmptyState()
    trajectory = Trajectory([state, ∅])
    @testset "Trajectory" begin
        @test trajectory == [state, ∅]
        @test trajectory[1] === state
    end
    metadata = Metadata(0.1)
    @testset "Metadata" begin
        @test metadata.weight === 0.1
    end
    particle = Particle(trajectory, metadata)
    @testset "Particle" begin
        @test particle.trajectory === trajectory
        @test particle.metadata === metadata
        particle.metadata.weight *= 2.0
        @test particle.metadata.weight === 0.2
    end
    params = Parameters(1.0, 0.1, 0.97, 0.03)
    @testset "Parameters" begin
        @test params.T === 1.0
        @test params.q === 0.1
        @test params.ps === 0.97
        @test params.pb === 0.03
    end
end

@testset "Motion Model" begin
    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    params = Parameters(1.0, 0.0, 0.97, 0.03)
    @test qk(state, params) == State(1, 2.0, 8.0, 10.0, 5.0, 6.0)
    params = Parameters(1.0, 0.1, 0.97, 0.03)
    @test logfk(state, state, params) === -0.0
end

end # testset
