module Obsea

using Interpolations

export State, EmptyState, Trajectory, Metadata, Particle, Parameters
export Grid, Scan, logl
export transition, move, birth, logf

include("particle.jl")
include("update.jl")
include("predict.jl")

end # module
