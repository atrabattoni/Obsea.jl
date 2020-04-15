using Obsea
using Test

@testset "State and Particle structures" begin
    state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    @test state.model === 1
    @test state.frequency === 2.0
    @test state.x === 3.0
    @test state.y === 4.0
    @test state.vx === 5.0
    @test state.vy === 6.0
    trajectory = [state]
    particle = Particle(trajectory, -1.0, 0.5)
    @test particle.trajectory === trajectory
    @test particle.logLikelihood === -1.0
    @test particle.weight === 0.5
end
