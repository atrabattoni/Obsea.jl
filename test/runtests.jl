using Obsea
using Test

@testset "Obsea.jl" begin

@testset "State, Trajectory, Metadata and Particle structures" begin

    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    @test state.model === 1
    @test state.frequency === 2.0
    @test state.x === 3.0
    @test state.y === 4.0
    @test state.vx === 5.0
    @test state.vy === 6.0

    trajectory = Trajectory([state, state])
    @test trajectory == [state, state]
    @test trajectory[1] === state

    metadata = Metadata(0.1)
    @test metadata.weight === 0.1

    particle = Particle(trajectory, metadata)
    @test particle.trajectory === trajectory
    @test particle.metadata === metadata
    particle.metadata.weight *= 2.0
    @test particle.metadata.weight === 0.2

end

@testset "Motion Model" begin
    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    @test qk(state, 0.0, 1.0) == State(1, 2.0, 8.0, 10.0, 5.0, 6.0)
end

end # testset
