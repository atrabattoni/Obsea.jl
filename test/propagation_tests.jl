import Obsea: tlim, toa, tdoa, TdoaLut

@testset "propagation" begin

    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    @unpack nmode, depth, celerity, ic, ib = propa

    @testset "tlim" begin
        @test tlim(0.0, 1, depth, celerity) == depth / celerity
        @test tlim(0.0, 2, depth, celerity) == 3 * depth / celerity
        @test tlim(π / 4, 1, depth, celerity) ≈ sqrt(2) * depth / celerity
    end

    @testset "toa" begin
        @test isnan(toa(0.0, 1, depth, celerity, ic, ib))
        @test isnan(toa(100_000.0, 1, depth, celerity, ic, ib))
        @test toa(depth, 1, depth, celerity, ic, ib) ==
              sqrt(2) * depth / celerity
        @test toa(depth, 2, depth / 3, celerity, ic, ib) ==
              sqrt(2) * depth / celerity
    end

    @testset "tdoa" begin
        @test isnan(tdoa(0.0, 2, depth, celerity, ic, ib))
        @test isnan(tdoa(100_000.0, 1, depth, celerity, ic, ib))
        @test tdoa(0.0, 1, depth, celerity, 0.0, ib) ≈ 2 * depth / celerity
    end

    @testset "TdoaLut" begin
        tdoalut = TdoaLut(propa, grid)
    end

end
