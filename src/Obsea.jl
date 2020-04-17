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
    init,
    Parameters
export Grid, Scan, update!, logl
export argsample, resample
export predict!, transition, move, birth, logf

include("types.jl")
include("predict.jl")
include("update.jl")
include("resample.jl")

end # module
