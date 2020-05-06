using Obsea
using Test

using Parameters
using StructArrays

@testset "Obsea.jl" begin
    include("utils_tests.jl")
    include("parameters_tests.jl")
    include("propagation_tests.jl")
    include("likelihood_tests.jl")
    include("particle_tests.jl")
    include("predict_tests.jl")
    include("update_tests.jl")
    include("resamplemove_tests.jl")
    include("track_tests.jl")
end
