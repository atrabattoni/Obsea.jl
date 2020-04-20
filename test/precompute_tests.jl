@testset "precompute" begin

    nmode = 3
    depth = 4340.0
    celerity = 1502.0

    @testset "tlim" begin
        import Obsea.tlim
        @test length(tlim(deg2rad(45.0), nmode, depth, celerity)) == nmode + 1
        @test tlim(deg2rad(0.0), nmode, depth, celerity) ==
              [(2i - 1) * depth / celerity for i = 1:nmode+1]
    end

    @testset "tdoa" begin
        import Obsea.tdoa
        tc = tlim(deg2rad(30.0), nmode, depth, celerity)
        tb = tlim(deg2rad(80.0), nmode, depth, celerity)
        @test all(isnan.(tdoa(0.0, nmode, depth, celerity, tc, tb)))
        @test all(isnan.(tdoa(90.0, nmode, depth, celerity, tc, tb)))
    end

    @testset "convsame" begin
        import Obsea.convsame
        @test convsame(ones(5), ones(3)) ≈ [2.0, 3.0, 3.0, 3.0, 2.0]
    end

end
