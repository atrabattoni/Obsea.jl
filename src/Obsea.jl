module Obsea

using DSP,
    Interpolations,
    LoopVectorization,
    Parameters,
    SpecialFunctions,
    StructArrays
import Pkg.TOML, Random.shuffle!

include("utils.jl")
include("parameters.jl")
include("propagation.jl")
include("likelihood.jl")
include("particle.jl")
include("predict.jl")
include("update.jl")
include("resamplemove.jl")
include("track.jl")

end # module
