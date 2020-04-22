@testset "azimuth" begin

    @testset "wrapcauchy" begin
        import Obsea.wrapcauchy
        @test_throws AssertionError wrapcauchy(0.0, 0.0, 0.0)
        @test_throws AssertionError wrapcauchy(0.0, 0.0, 1.0)
        @test wrapcauchy(0.0, 0.0, 1 / sqrt(2)) ≈ 0.5 / (1 - 1 / sqrt(2))^2
    end

    @testset "rollprod" begin
        import Obsea.rollprod
        x = fill(2.0, (3, 2))
        n = 3
        @test rollprod(x, n) == [
            4.0 4.0
            8.0 8.0
            4.0 4.0
        ]
    end

    @testset "Azimuth" begin
        import Obsea.Azimuth
        @test_nowarn Azimuth(0.9, 5, 0.0:3.0:360)
    end

    @testset "precomp" begin
        import Obsea.precomp
        z = zeros(200)
        arange = range(0, 360, length=121)
        model= Azimuth(0.9, 5, arange)
        @test size(precomp(z, model)) == (200, 121, 2)
    end

end
