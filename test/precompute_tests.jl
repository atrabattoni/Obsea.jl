import Obsea: precompute, distribution

@testset "precompute" begin

      dict = TOML.parsefile("parameters.toml")
      models, propa, grid = parameters(dict, 50.0, 1024)
      tdoalut = TDOALUT(propa, grid)
      zr = zeros(grid.Nτ)
      lams = [models[i].lam for i = 1:grid.Nm]
      @assert all(x -> x == lams[1], lams)
      lam = first(lams)
      Nmode = propa.Nmode

      @unpack Nmode, depth, celerity, ic, ib, sigma = propa
      propa = Propagation(
            Nmode = Nmode,
            depth = depth,
            celerity = celerity,
            ic = 0.0,
            ib = 90.0,
            sigma = sigma,
      )
      tdoalut = TDOALUT(propa, grid)
      ℓ = precompute(zr, tdoalut, models, grid)
      @test ℓ ≈ fill(exp(-lam / 2)^Nmode, grid.Nr, grid.Nm)

      @unpack Nmode, depth, celerity, ic, ib, sigma = propa
      propa = Propagation(
            Nmode = Nmode,
            depth = depth,
            celerity = celerity,
            ic = 90.0,
            ib = 0.0,
            sigma = sigma,
      )
      tdoalut = TDOALUT(propa, grid)
      ℓ = precompute(zr, tdoalut, models, grid)
      @test ℓ ≈ fill(1.0, grid.Nr, length(1:grid.Nm))

      za = zeros(grid.Nf)
      ℓ = precompute(za, models, grid)
      @test length(ℓ) == grid.Nf * grid.Na * length(1:grid.Nm)

      @test_nowarn precompute(zr, za, tdoalut, models, grid)

end
