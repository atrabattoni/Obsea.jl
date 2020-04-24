import Obsea: resample

@testset "resample" begin

    @testset "resample" begin
        particle_1 = [State(1, 0, 0, 0, 0, 0)]
        particle_2 = [State(2, 0, 0, 0, 0, 0)]
        weights = [0.0, 1.0]
        cloud = [particle_1, particle_2]
        weights, cloud = resample(weights, cloud)
        @test cloud[1] == particle_2
        @test cloud[2] == particle_2
        @test weights[1] == 0.5
        @test weights[2] == 0.5
    end

end
