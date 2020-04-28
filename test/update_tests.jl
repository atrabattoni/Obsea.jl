import Obsea: likelihoodratio, update!

@testset "update" begin

    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    tdoalut = TDOALUT(propa, grid)
    zr = zeros(grid.Nτ)
    za = zeros(grid.Nf)
    ℓ = precompute(zr, za, tdoalut, models, grid)
    ∅ = State()
    ship1 = State(1, 5.0, 1000.0, 0.0, 5.0, 5.0)
    ship2 = State(1, 15.0, 1000.0, 1000.0, 5.0, 5.0)
    whale = State(2, 15.0, 1000.0, 0.0, 5.0, 5.0)

    @testset "likelihoodratio" begin
        @test likelihoodratio(ℓ, ∅, models, grid) === 1.0
        @test abs(likelihoodratio(ℓ, ship1, models, grid) - 1) < 10 * eps(1.0)
        @test likelihoodratio(ℓ, ship2, models, grid) < 1.0
        @test likelihoodratio(ℓ, whale, models, grid) < 1.0
    end

    @testset "update" begin
        import Obsea.update!
        weights = [1.0]
        cloud = [[whale]]
        update!(weights, cloud, ℓ, models, grid)
        @test weights[1] ≈ 1.0

    end

end
