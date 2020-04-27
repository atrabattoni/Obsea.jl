import Obsea: precompute

@testset "precompute" begin

      dict = TOML.parsefile("parameters.toml")
      models, propa, grid = parameters(dict, 50.0, 1024)
      tdoalut = TDOALUT(propa, grid)
      zr = zeros(grid.Nτ)
      lams = [models[i].lam for i in 1:grid.Nm]
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
      ℓ, itp = precompute(zr, tdoalut, models, grid)
      @test ℓ ≈
            fill(exp(-lam / 2)^Nmode, grid.Nr, grid.Nm)
      @test itp.(collect(grid.r), fill(1, grid.Nr)) ≈ ℓ[:, 1]
      @test itp.(collect(grid.r), fill(2, grid.Nr)) ≈ ℓ[:, 2]

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
      ℓ, itp = precompute(zr, tdoalut, models, grid)
      @test ℓ ≈ fill(1.0, grid.Nr, length(1:grid.Nm))
      @test itp.(collect(grid.r), fill(1, grid.Nr)) ≈ ℓ[:, 1]
      @test itp.(collect(grid.r), fill(2, grid.Nr)) ≈ ℓ[:, 2]

      za = zeros(grid.Nf)
      ℓ, itp = precompute(za, models, grid)
      @test length(ℓ) ==
            grid.Nf * grid.Na * length(1:grid.Nm)
      @test itp(rand(grid.f), rand(grid.a), rand(1:grid.Nm)) > 0

      @test_nowarn precompute(zr, za, tdoalut, models, grid)

end
