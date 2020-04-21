@testset "range" begin

    nmode = 3
    depth = 4340.0
    celerity = 1502.0
    ic = deg2rad(10.0)
    ib = deg2rad(80.0)

    lam = 1.0
    rrange = range(0.0, 30_000.0, step = 100.0)
    τrange = range(0.0, 10.24, length = 513)
    z = zeros(513)

    @testset "tlim" begin
        import Obsea.tlim
        @test tlim(0.0, 1, depth, celerity) == depth / celerity
        @test tlim(0.0, 2, depth, celerity) == 3 * depth / celerity
        @test tlim(π / 4, 1, depth, celerity) ≈ sqrt(2) * depth / celerity
    end

    @testset "toa" begin
        import Obsea.toa
        @test isnan(toa(0.0, 1, depth, celerity, ic, ib))
        @test isnan(toa(100_000.0, 1, depth, celerity, ic, ib))
        @test toa(depth, 1, depth, celerity, ic, ib) ==
              sqrt(2) * depth / celerity
        @test toa(depth, 2, depth / 3, celerity, ic, ib) ==
              sqrt(2) * depth / celerity
    end

    @testset "tdoa" begin
        import Obsea.tdoa
        @test isnan(tdoa(0.0, 2, depth, celerity, ic, ib))
        @test isnan(tdoa(100_000.0, 1, depth, celerity, ic, ib))
        @test tdoa(0.0, 1, depth, celerity, 0.0, ib) ≈ 2 * depth / celerity
    end

    @testset "propagation" begin
        import Obsea.propagation
        @test propagation([0.0], 2, depth, celerity, 0.0, π / 2) ≈
              [2 * depth / celerity 2 * depth / celerity]
    end

    @testset "window" begin
        import Obsea.window
        @test sum(window(7, 0.25)) ≈ 1.0
        @test window([5, 7], [0.5, 0.25])[1] == window(5, 0.5)
    end

    @testset "convsame" begin
        import Obsea.convsame
        u = rand(10)
        @test convsame(u, ones(1)) == u
        @test convsame(ones(5), ones(3)) ≈ [2.0, 3.0, 3.0, 3.0, 2.0]
    end

    @testset "Range" begin
        import Obsea.Range
        n = fill(1, nmode)
        σ = fill(1.0, nmode)
        model = Range(nmode, depth, celerity, ic, ib, lam, n, σ, rrange, τrange)
        @test length(model.v) == nmode
        @test size(model.τ) == (length(rrange), nmode)
    end

    @testset "precomp" begin
        import Obsea.precomp
        n = fill(1, nmode)
        σ = fill(1.0, nmode)
        model =
            Range(nmode, depth, celerity, 0, pi / 2, lam, n, σ, rrange, τrange)
        @test precomp(z, model) ≈ fill(exp(-lam / 2)^nmode, length(rrange))
        model =
            Range(nmode, depth, celerity, pi / 2, 0, lam, n, σ, rrange, τrange)
        @test precomp(z, model) ≈ fill(1.0, length(rrange))
    end
end
