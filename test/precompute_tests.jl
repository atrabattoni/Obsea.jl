@testset "precompute" begin

    nmode = 3
    depth = 4340.0
    celerity = 1502.0

    @testset "tlim" begin
        import Obsea.tlim
        @test length(tlim(deg2rad(45.0), depth, celerity, nmode)) == nmode + 1
        @test tlim(deg2rad(0.0), depth, celerity, nmode) ==
              [(2i - 1) * depth / celerity for i=1:nmode + 1]
    end

    @testset "tdoa" begin
        import Obsea.tdoa
        tc = tlim(deg2rad(30.0), depth, celerity, nmode)
        tb = tlim(deg2rad(80.0), depth, celerity, nmode)
        @test all(isnan.(tdoa(0.0, depth, celerity, nmode, tc, tb)))
        @test all(isnan.(tdoa(90.0, depth, celerity, nmode, tc, tb)))
    end

end
