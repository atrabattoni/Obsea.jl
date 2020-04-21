@testset "update.jl" begin

    grid = Grid(
        (0.0:100.0:30000),
        range(0.0, 25.0, length = 513),
        range(0.0, 2π, length = 121),
        (1:2),
    )
    scan = Scan(
        ones(length(grid.range_r)),
        ones(length(grid.range_f), length(grid.range_a), 2),
        grid,
    )
    ∅ = EmptyState()
    ship = State(1, 5.0, 1000.0, 1000.0, 0.0, 0.0)
    whale = State(2, 5.0, 1000.0, 1000.0, 0.0, 0.0)
    params = Parameters(1.0, 0.0, 0.97, 0.03, 0.5, grid)

    @testset "Grid" begin
        @test convert(Float64, grid.range_r.step) == 100.0
        @test grid.range_f.len === 513
        @test convert(Float64, grid.range_f.step) == 25.0 / 512
        @test grid.range_a[end] === 2π
        @test collect(grid.range_m) == [1, 2]
    end

    @testset "Scan" begin
        @test scan.itp_r(1100.0) ≈ 1.0
        @test scan.itp_fam(10.0, π / 3, 1) ≈ 1.0
    end

    @testset "logl" begin
        import Obsea.logl
        @test logl(scan, ∅, params) == 0.0
        @test abs(logl(scan, ship, params)) < 1e-15
        @test abs(logl(scan, whale, params)) < 1e-15
    end

    @testset "update" begin
        import Obsea.update!
        weights = [1.0]
        cloud = [[ship]]
        update!(weights, cloud, scan, params)
        @test weights[1] ≈ 1.0

    end

end
