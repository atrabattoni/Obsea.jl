import Obsea: LikelihoodSlice, likelihood, marginalize

@testset "likelihood, marginalize" begin

      models, propa, grid = parameters("parameters.toml", 50.0, 1024)
      tdoalut = TDOALUT(propa, grid)
      Nt = 10
      zr = zeros(grid.Nτ, Nt)
      lams = models.lam
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
      ℓr = likelihood(zr, tdoalut, models, grid)
      @test ℓr ≈ fill(exp(-lam / 2)^Nmode, grid.Nr, Nt, grid.Nm)

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
      ℓ = likelihood(zr, tdoalut, models, grid)
      @test ℓ ≈ fill(1.0, grid.Nr, Nt, length(1:grid.Nm))

      za = zeros(grid.Nf, Nt)
      ℓa = likelihood(za, models, grid)
      @test length(ℓa) == grid.Nf * grid.Na * Nt *length(1:grid.Nm)

      ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)

end
