using Obsea
using Test

using Parameters
using Pkg.TOML

@testset "Obsea.jl" begin
    @time include("utils_tests.jl")
    @time include("parameters_tests.jl")
    @time include("propagation_tests.jl")
    @time include("precompute_tests.jl")
    @time include("particle_tests.jl")
    @time include("predict_tests.jl")
    # include("update_tests.jl")
    # include("reamplemove_tests.jl")
    # include("tracker_tests.jl")
end
