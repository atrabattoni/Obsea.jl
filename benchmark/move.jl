using BenchmarkTools
using Profile
using Parameters
using LoopVectorization
import Obsea: State

function kinematic!(newcloud, cloud, T, Q)
    dims = size(cloud)
    q = Q[cloud.m]
    ax = q .* randn(dims...)
    ay = q .* randn(dims...)
    newcloud.x .= cloud.x .+ cloud.vx .* T .+ ax .* (T^2 / 2.0)
    newcloud.y .= cloud.x .+ cloud.vy .* T .+ ay .* (T^2 / 2.0)
    newcloud.vx .= cloud.vx .+ ax .* T
    newcloud.vy .= cloud.vy .+ ay .* T
    return newcloud
end

function kill!(cloud, Ps)
    dims = size(cloud)
    ps = Ps[cloud.m]
    mask = rand(dims...) .>= ps
    @views fill!(cloud[mask], State())
    return cloud
end


function move!(newcloud, cloud, T, Ps, Q)
    kinematic!(newcloud, cloud, T, Q)
    kill!(newcloud, Ps)
    return newcloud
end


function wrapmove!(cloud, mask, T, Ps, Q)
    cloud[mask] .= move!(cloud[mask], T, Ps, Q)
end


Nm = 2
N = 1000
Ps = [0.9, 0.8]
Q = [0.1, 0.2]
T = 10.0
newcloud = StructVector{State}((
    rand(1:Nm, N),
    rand(N),
    rand(N),
    rand(N),
    rand(N),
    rand(N),
))
cloud = deepcopy(cloud)
mask = rand(N) .< 0.5

println()
println("reference")
@btime move!(newcloud, $cloud, $T, $Ps, $Q) setup =
    (newcloud = deepcopy($newcloud)) evals = 1
println("view")
@btime move!(view(c, $mask), T, Ps, Q) setup = (c = deepcopy($cloud)) evals = 1
println("wrap")
@btime wrapmove!(c, $mask, $T, $Ps, $Q) setup = (c = deepcopy($cloud)) evals = 1
println()
c = deepcopy(cloud)
wrapmove!(c, mask, T, Ps, Q)

# @profiler for i = 1:1000 #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Bottleneck = randn
#     c = deepcopy(cloud)
#     move!(c, cloud, T, Ps, Q)
# end
