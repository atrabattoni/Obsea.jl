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
    ship = ShipState(5.0, 1000.0, 1000.0, 0.0, 0.0)
    whale = WhaleState(5.0, 1000.0, 1000.0, 0.0, 0.0)
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
        @test logl(scan, ∅) == 0.0
        @test abs(logl(scan, ship, params)) < 1e-15
        @test abs(logl(scan, whale, params)) < 1e-15
    end

    @testset "update" begin
        metadata = Metadata(1.0)
        trajectory = Trajectory([ship])
        particle = Particle(trajectory, metadata)
        cloud = Cloud([particle])
        update!(cloud, scan, params)
        @test particle.metadata.weight ≈ 1.0

    end

end
