module Obsea

using Interpolations

export State,
    EmptyState,
    ShipState,
    WhaleState,
    Trajectory,
    Metadata,
    Particle,
    Cloud,
    Parameters
export Grid, Scan, update!, logl
export argsample, resample
export predict!, transition, move, birth, logf

include("particle.jl")
include("update.jl")
include("resample.jl")
include("predict.jl")

end # module
