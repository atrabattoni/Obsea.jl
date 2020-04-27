import Obsea.track

@testset "tracker" begin
    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
    nt = 10
    ceps = zeros(grid.Nτ, nt)
    az = zeros(grid.Nf, nt)
    track(ceps, az, "parameters.toml", 50.0, 1024)
end
