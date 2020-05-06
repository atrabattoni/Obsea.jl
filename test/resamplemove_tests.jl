import Obsea: resample!, alivelength

@testset "resample!" begin
    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    s0 = State(0, 0, 0, 0, 0, 0)
    s1 = State(1, 0, 0, 0, 0, 0)
    s2 = State(2, 0, 0, 0, 0, 0)
    weights = [2 / 3, 1 / 3, 0.0]
    particles = StructArray([
        s1 s1 s2
        s0 s2 s0
    ])
    ℓ = nothing
    resample!(weights, particles, ℓ, models, grid)
    @test particles[:, 1] == StructVector([s1, s0])
    @test particles[:, 2] == StructVector([s1, s2])
    @test particles[:, 3] == StructVector([s1, s0])
    @test weights == ones(3) ./ 3

end

@testset "alivelength" begin
    ∅ = State(0, 0, 0, 0, 0, 0)
    s = State(1, 0, 0, 0, 0, 0)
    @test alivelength(StructVector([∅, ∅, ∅, ∅, s, s, s])) == 3
    @test alivelength(StructVector([∅, ∅, ∅, ∅, s, s, ∅])) == 0
    @test alivelength(StructVector([∅, ∅, ∅, ∅, ∅, ∅, s])) == 1
    @test alivelength(StructVector([∅, ∅, ∅, ∅, ∅, ∅, ∅])) == 0
    @test alivelength(StructVector([s, s, s, s, s, s, s])) == 7
end
