using Obsea
using Test

@testset "Obsea.jl" begin
    include("particle_tests.jl")
    include("predict_tests.jl")
    include("update_tests.jl")
end
