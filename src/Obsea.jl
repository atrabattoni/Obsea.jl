module Obsea

using Interpolations

export State, EmptyState, Trajectory, Metadata, Particle, Parameters
export move, logfk

include("types.jl")
include("move.jl")
include("logfk.jl")

end # module
