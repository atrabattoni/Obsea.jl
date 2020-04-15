using Obsea
using Test

@testset "State" begin
    s = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
    @test s.m === 1
    @test s.x === 2.0
    @test s.y === 3.0
    @test s.vx === 4.0
    @test s.vy === 5.0
end
