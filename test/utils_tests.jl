import Obsea:
    window, limit, symbolize, convsame, wrapcauchy, rollprod, argsample, ra2xy

@testset "utils" begin

    @testset "window" begin
        @test window(1, 0.5) ≈ [1.0]
        @test sum(window(7, 0.25)) ≈ 1.0
    end

    @testset "limit" begin
        @test_throws AssertionError limit(4.0:-1.0, 1.0, 3.0)
        @test_throws AssertionError limit(-1.0:4.0, 3.0, 1.0)
        @test limit(-1.0:4.0, 1.0, 3.0) == 1.0:3.0
        @test limit(-1.0:4.0, 0.5, 3.5) == 1.0:3.0
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
    end

    @testset "wrapcauchy" begin
        @test_throws AssertionError wrapcauchy(0.0, 0.0, 0.0)
        @test_throws AssertionError wrapcauchy(0.0, 0.0, 1.0)
        @test wrapcauchy(0.0, 0.0, 1 / sqrt(2)) ≈ 0.5 / (1 - 1 / sqrt(2))^2
    end

    @testset "rollprod" begin
        x = fill(2.0, (3, 2))
        @test_throws AssertionError rollprod(x, 2)
        @test rollprod(x, 3) == [
            4.0 4.0
            8.0 8.0
            4.0 4.0
        ]
    end

    @testset "argsample" begin
        @test argsample([0.0, 1.0]) === 2
        @test argsample([0.0, 0.0]) === 0
        @test argsample(ones(10) / 10, 10) == collect(1:10)
        @test argsample([0.0, 0.0], 10) == fill(0, 10)
    end

end
