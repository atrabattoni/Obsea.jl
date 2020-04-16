module Obsea

using Interpolations

export State, EmptyState, Trajectory, Metadata, Particle, Parameters
export Grid, Scan, logl
export argsample
export transition, move, birth, logf

include("particle.jl")
include("update.jl")
include("resample.jl")
include("predict.jl")

end # module
