import Obsea: make, update!

@testset "update" begin

    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    tdoalut = TDOALUT(propa, grid)
    zr = zeros(grid.Nτ, 1)
    za = zeros(grid.Nf, 1)
    ℓr = likelihood(zr, tdoalut, models, grid)
    ℓa = likelihood(za, models, grid)
    ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)
    @views ℓ = LikelihoodSlice(ℓr[:, 1, :], ℓa[:, :, 1, :], ℓm[1, :], ℓΣm[1])
    ∅ = State()
    ship1 = State(1, 5.0, 1000.0, 0.0, 5.0, 5.0)
    ship2 = State(1, 15.0, 1000.0, 1000.0, 5.0, 5.0)
    whale = State(2, 15.0, 1000.0, 0.0, 5.0, 5.0)


    @testset "update" begin
        weights = [1.0]
        cloud = StructVector([whale])
        update!(weights, cloud, ℓ, models, grid)
        @test weights[1] ≈ 1.0
    end

end
