import Obsea: resample!

@testset "resample!" begin
    s0 = State(0, 0, 0, 0, 0, 0)
    s1 = State(1, 0, 0, 0, 0, 0)
    s2 = State(2, 0, 0, 0, 0, 0)
    weights = [2 / 3, 1 / 3, 0.0]
    particles = StructArray([
        s1 s1 s2
        s0 s2 s0
    ])
    resample!(weights, particles)
    @test particles[:, 1] == StructVector([s1, s0])
    @test particles[:, 2] == StructVector([s1, s2])
    @test particles[:, 3] == StructVector([s1, s0])
    @test weights == ones(3) ./ 3

end
