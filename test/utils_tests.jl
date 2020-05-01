import Obsea:
    window,
    limit,
    grid2range,
    symbolize,
    convsame!,
    convsame,
    rollprod!,
    rollprod,
    argsample,
    ra2xy,
    xy2ra,
    counts

@testset "utils" begin

    @testset "window" begin
        @test window(0.1) ≈ [1.0]
        σ = 2.5
        w = window(σ)
        Nv = length(w)
        @test isodd(Nv)
        @test sum(w) ≈ 1.0
        x = collect(-(Nv ÷ 2):(Nv÷2))
        @test sum(x .* w) ≈ 0.0
        @test σ^2 / 2 < sum((x .^ 2) .* w) / sum(w) < σ^2
    end

    @testset "limit" begin
        @test_throws AssertionError limit(4.0:-1.0, 1.0, 3.0)
        @test_throws AssertionError limit(-1.0:4.0, 3.0, 1.0)
        @test limit(-1.0:4.0, 1.0, 3.0) == 1.0:3.0
        @test limit(-1.0:4.0, 0.5, 3.5) == 1.0:3.0
    end

    @testset "grid2range" begin
        @test grid2range(0.2:0.1:3) == 0.2:0.1:3
    end

    @testset "symbolize" begin
        @test symbolize(Dict("x" => 1, "y" => 2)) == Dict(:x => 1, :y => 2)
    end

    @testset "convsame" begin
        @test_throws AssertionError convsame(ones(3), ones(5))
        @test_throws AssertionError convsame(ones(5), ones(2))
        u = rand(10)
        @test convsame(u, ones(1)) == u
        @test convsame(ones(5), ones(3)) ≈ [2.0, 3.0, 3.0, 3.0, 2.0]
        out = similar(ones(5))
        convsame!(out, ones(5), ones(3))
        @test out ≈ [2.0, 3.0, 3.0, 3.0, 2.0]
    end

    @testset "rollprod" begin
        x = fill(2.0, 3)
        @test_throws AssertionError rollprod(x, 2)
        @test rollprod(x, 3) == [4.0, 8.0, 4.0]
        out = similar(x)
        rollprod!(out, x, 3)
        @test out == [4.0, 8.0, 4.0]
    end

    @testset "argsample" begin
        @test argsample([0.0, 1.0], 1, scale = 1) == [2]
        @test argsample([0.0, 0.0], 1, scale = 1) == [0]
        @test sort(argsample(ones(10) / 10, 10)) == collect(1:10)
        @test argsample([0.0, 0.0], 10, scale = 1.0) == fill(0, 10)
    end

    @testset "polar coordinate" begin
        x = randn()
        y = randn()
        @test all(ra2xy(xy2ra(x, y)...) .≈ (x, y))
    end

    @testset "counts" begin
        @test counts([1, 1, 2, 0, 1], 3) == Dict(0 => 1, 1 => 3, 2 => 1, 3 => 0)
    end

end
