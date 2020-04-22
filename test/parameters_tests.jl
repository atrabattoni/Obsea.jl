import Obsea: Propagation, Model, Grid, parameters

@testset "parameters" begin
    models, propa, grid = parameters("parameters.toml", 50.0, 1024)
end
