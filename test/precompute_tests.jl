import Obsea: CDF, ITP, precompute

@testset "precompute" begin

      dict = TOML.parsefile("parameters.toml")
      models, propa, grid = parameters(dict, 50.0, 1024)
      tdoalut = TDOALUT(propa, grid)
      zr = zeros(length(grid.τrange))
      lams = [models[i].lam for i in grid.mrange]
      @assert all(x -> x == lams[1], lams)
      lam = first(lams)
      nmode = propa.nmode

      @unpack nmode, depth, celerity, ic, ib, sigma = propa
      propa = Propagation(
            nmode = nmode,
            depth = depth,
            celerity = celerity,
            ic = 0.0,
            ib = 90.0,
            sigma = sigma,
      )
      tdoalut = TDOALUT(propa, grid)
      cdf, itp = precompute(zr, tdoalut, models, grid)
      y = fill(exp(-lam / 2)^nmode, length(grid.rrange), length(grid.mrange))
      @test itp.(collect(grid.rrange), fill(1, length(grid.rrange))) ≈ y[:, 1]
      @test itp.(collect(grid.rrange), fill(2, length(grid.rrange))) ≈ y[:, 2]
      @test cdf ≈ cumsum(vec(y)) / sum(y)

      @unpack nmode, depth, celerity, ic, ib, sigma = propa
      propa = Propagation(
            nmode = nmode,
            depth = depth,
            celerity = celerity,
            ic = 90.0,
            ib = 0.0,
            sigma = sigma,
      )
      tdoalut = TDOALUT(propa, grid)
      cdf, itp = precompute(zr, tdoalut, models, grid)
      y = fill(1.0, length(grid.rrange), length(grid.mrange))
      @test itp.(collect(grid.rrange), fill(1, length(grid.rrange))) ≈ y[:, 1]
      @test itp.(collect(grid.rrange), fill(2, length(grid.rrange))) ≈ y[:, 2]
      @test cdf ≈ cumsum(vec(y)) / sum(y)

      za = zeros(length(grid.frange))
      cdf, itp = precompute(za, models, grid)
      @test length(cdf) ==
            length(grid.frange) * length(grid.arange) * length(grid.mrange)
      @test itp(rand(grid.frange), rand(grid.arange), rand(grid.mrange)) > 0

      @test_nowarn precompute(zr, za, tdoalut, models, grid)

end
