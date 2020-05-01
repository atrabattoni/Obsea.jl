import Obsea: Propagation, Model, Grid, parameters

@testset "parameters" begin
    models, propa, grid = parameters("parameters.toml", 50.0, 1024)

    @test models.q == [0.1, 0.1]
    @test models.vmin == [3.0, 0.0]
    @test models.vmax == [12.0, 6.0]
    @test models.ps == [0.99, 0.99]
    @test models.pd == [1.0, 0.2]
    @test models.lam == [1.0, 1.0]
    @test models.mrl == [0.9, 0.9]
    @test models.n == [1, 5]

    @test propa.Nmode == 3
    @test propa.depth == 4340.0
    @test propa.celerity == 1502.0
    @test propa.ic == deg2rad(10.0)
    @test propa.ib == deg2rad(80.0)
    @test propa.sigma == [1, 3, 5]

    @test grid.r == collect(0.0:100.0:30000.0)
    @test first(grid.f) > 11.0
    @test last(grid.f) < 24.0
    @test length(grid.a) == 121
    @test first(grid.a) == 0.0
    @test last(grid.a) == 2π
    @test length(grid.τ) == 513
    @test first(grid.τ) == 0.0
    @test last(grid.τ) == 10.24

end
