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

    @testset "precomp_r" begin
        import Obsea.precomp_r
        z = rand(513).^2
        lam = 1.0
        nmode = 3
        range_r = 0.0:100.0:30000.0
        range_tau = range(0.0, 10.24, length = 513)
        n = [1, 3, 5]
        σ = 0.5 / 2.0
        tc = tlim(deg2rad(30.0), nmode, depth, celerity)
        tb = tlim(deg2rad(80.0), nmode, depth, celerity)
        precomp_r(
            z,
            lam,
            nmode,
            range_r,
            range_tau,
            n,
            σ,
            depth,
            celerity,
            tc,
            tb,
        )
    end
end
