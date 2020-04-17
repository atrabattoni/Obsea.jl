@testset "resample.jl" begin

    particle_1 = Particle(Trajectory([ShipState(0, 0, 0, 0, 0)]), Metadata(0.0))
    particle_2 = Particle(Trajectory([WhaleState(0, 0, 0, 0, 0)]), Metadata(1.0))
    cloud = Cloud([particle_1, particle_2])

    @testset "argsample" begin
        w = ones(10)
        w /= sum(w)
        @test argsample(w, 10) == collect(1:10)
    end

    @testset "resample" begin
        cloud = resample(cloud)
        @test typeof(cloud[1].trajectory[1]) == WhaleState
        @test typeof(cloud[2].trajectory[1]) == WhaleState
        @test cloud[1].metadata.weight == 0.5
        @test cloud[2].metadata.weight == 0.5
    end

end
