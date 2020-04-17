using Obsea
using Test

@testset "Obsea.jl" begin
    include("types_tests.jl")
    include("predict_tests.jl")
    include("update_tests.jl")
    include("resample_tests.jl")
end
