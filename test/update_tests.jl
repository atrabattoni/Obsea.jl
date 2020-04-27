import Obsea: likelihoodratio, update!

@testset "update" begin

    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    tdoalut = TDOALUT(propa, grid)
    zr = zeros(grid.Nτ)
    za = zeros(grid.Nf)
    ℓ, itp = precompute(zr, za, tdoalut, models, grid)
    ∅ = State()
    ship1 = State(1, 5.0, 1000.0, 0.0, 5.0, 5.0)
    ship2 = State(1, 15.0, 1000.0, 1000.0, 5.0, 5.0)
    whale = State(2, 15.0, 1000.0, 0.0, 5.0, 5.0)

    @testset "likelihoodratio" begin
        @test likelihoodratio(itp, ∅, models) === 1.0
        @test abs(likelihoodratio(itp, ship1, models) - 1) < 10 * eps(1.0)
        @test likelihoodratio(itp, ship2, models) < 1.0
        @test likelihoodratio(itp, whale, models) < 1.0
    end

    @testset "update" begin
        import Obsea.update!
        weights = [1.0]
        cloud = [[whale]]
        update!(weights, cloud, itp, models)
        @test weights[1] ≈ 1.0

    end

end
