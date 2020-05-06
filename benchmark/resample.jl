using BenchmarkTools
using Profile
using Parameters
using LoopVectorization
using StructArrays
import Obsea: State, argsample

function resample(weights, particles, idx)
    @views particles = deepcopy(particles[:, idx])
    weights = fill(1 / length(weights), length(weights))
    return weights, particles
end

function resample!(weights, particles, idx)
    Np = length(weights)
    keep = Vector{Bool}(undef, Np)
    stack = Vector{Int}(undef, 0)
    for j in idx
        if !keep[j]
            keep[j] = true
        else
            push!(stack, j)
        end
    end
    for j = 1:Np
        if !keep[j]
            k = pop!(stack)
            @views particles[:, j] .= particles[:, k]
        end
    end
    fill!(weights, 1 / Np)
    return weights, particles
end

function resample2!(keep, stack, weights, particles, idx)
    Np = length(weights)
    count = 0
    for j in idx
        if !keep[j]
            keep[j] = true
        else
            count += 1
            stack[count] = j
        end
    end
    for j = 1:Np
        if !keep[j]
            k = stack[count]
            count -= 1
            replpacefieldscolumn!(particles, j, k)
        end
    end
    #fill!(weights, 1 / Np)
    return weights, particles
end

function replpacefieldscolumn!(particles, i, j)
    replacecolumn!(particles.m, i, j)
    replacecolumn!(particles.f, i, j)
    replacecolumn!(particles.x, i, j)
    replacecolumn!(particles.y, i, j)
    replacecolumn!(particles.vx, i, j)
    replacecolumn!(particles.vy, i, j)

end

function replacecolumn!(A, i, j)
    @inbounds for k = 1:size(A, 1)
        A[k, i] = A[k, j]
    end
end

Nm = 2
Np = 5000
Nt = 5000
weights = rand(Np) .^ 2
weights /= sum(weights)
particles = StructArray{State}((
    rand(1:Nm, Nt, Np),
    rand(Nt, Np),
    rand(Nt, Np),
    rand(Nt, Np),
    rand(Nt, Np),
    rand(Nt, Np),
))
idx = argsample(weights, length(weights); scale = 1.0)


println()
println("reference")
@btime resample($weights, $particles, $idx)
println("in-place")
@btime resample!(w, p, $idx) setup =
    (w = copy($weights); p = deepcopy($particles)) evals = 1
println("in-total-place")
@btime resample2!(k, s, w, p, $idx) setup = (
    k = Vector{Bool}(undef, Np); s = Vector{Int}(undef, Np); w = copy($weights); p = deepcopy($particles)
) evals = 1
# @code_warntype resample!(copy(weights), particles, idx)
# @profiler for i = 1:100
#     resample!(copy(weights), particles, idx)
# end
