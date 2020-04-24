import Obsea: precompute

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
      ℓ, itp = precompute(zr, tdoalut, models, grid)
      @test ℓ ≈
            fill(exp(-lam / 2)^nmode, length(grid.rrange), length(grid.mrange))
      @test itp.(collect(grid.rrange), fill(1, length(grid.rrange))) ≈ ℓ[:, 1]
      @test itp.(collect(grid.rrange), fill(2, length(grid.rrange))) ≈ ℓ[:, 2]

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
      ℓ, itp = precompute(zr, tdoalut, models, grid)
      @test ℓ ≈ fill(1.0, length(grid.rrange), length(grid.mrange))
      @test itp.(collect(grid.rrange), fill(1, length(grid.rrange))) ≈ ℓ[:, 1]
      @test itp.(collect(grid.rrange), fill(2, length(grid.rrange))) ≈ ℓ[:, 2]

      za = zeros(length(grid.frange))
      ℓ, itp = precompute(za, models, grid)
      @test length(ℓ) ==
            length(grid.frange) * length(grid.arange) * length(grid.mrange)
      @test itp(rand(grid.frange), rand(grid.arange), rand(grid.mrange)) > 0

      @test_nowarn precompute(zr, za, tdoalut, models, grid)

end
