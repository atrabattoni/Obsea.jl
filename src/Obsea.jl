module Obsea

using DSP
using Interpolations
using Parameters
using SpecialFunctions

export State,
    EmptyState,
    ShipState,
    WhaleState,
    Trajectory,
    Metadata,
    Particle,
    Cloud,
    init,
    Parameters
export Grid, Scan, update!, logl
export argsample, resample
export predict!, transition, move, birth, logf

include("types.jl")
include("range.jl")
include("predict.jl")
include("update.jl")
include("resample.jl")

end # module
