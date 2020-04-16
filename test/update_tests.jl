@testset "update.jl" begin

    @testset "Grid" begin
        grid = Grid(
            (0.0:1000.0:30000),
            range(0.0, 25.0, length = 513),
            (0.0:5.0:360),
            (1:2),
        )
        @test convert(Float64, grid.range_r.step) == 1000.0
        @test grid.range_f.len === 513
        @test convert(Float64, grid.range_f.step) == 25.0 / 512
        @test grid.range_a[end] === 360.0
        @test collect(grid.range_m) == [1, 2]
    end

    @testset "Scan" begin
        grid = Grid(
            range(0.0, 1000.0, length = 5),
            range(0.0, 10.0, length = 7),
            range(0.0, 360.0, length = 9),
            (1:2),
        )
        scan = Scan(ones(5), ones(7, 9, 2), grid)
        @test scan.itp_r(500.0) ≈ 1.0
        @test scan.itp_fam(5.0, 120.0, 1) ≈ 1.0
    end

    @testset "logl" begin
        grid = Grid(
            range(0.0, 1000.0, length = 5),
            range(0.0, 10.0, length = 7),
            range(0.0, 2π, length = 9),
            (1:2),
        )
        scan = Scan(ones(5), ones(7, 9, 2), grid)
        ∅ = EmptyState()
        params = Parameters(1.0, 0.0, 0.97, 0.03, 0.5)
        @test logl(scan, ∅) == 0.0
        state = State(1, 5.0, 1000.0, 1000.0, 0.0, 0.0)
        @test abs(logl(scan, state, params)) < 1e-15
        state = State(2, 5.0, 1000.0, 1000.0, 0.0, 0.0)
        @test abs(logl(scan, state, params)) < 1e-15
    end

end
