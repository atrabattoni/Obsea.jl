import Obsea: track, fetch

@testset "tracker" begin
    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    Nt = 10
    zr = zeros(grid.Nτ, Nt)
    za = zeros(grid.Nf, Nt)
    track(zr, za, "parameters.toml", 50.0, 1024, 100)
end
