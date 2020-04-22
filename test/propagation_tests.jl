import Obsea: tlim, toa, tdoa

@testset "propagation" begin

    nmode = 3
    depth = 4340.0
    celerity = 1502.0
    ic = deg2rad(10.0)
    ib = deg2rad(80.0)
    
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

end
