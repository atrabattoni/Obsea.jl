function resample!(weights, particles)
    idx = argsample(weights, length(weights); scale = 1.0)
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
            @views particles[:, j] .= particles[:, pop!(stack)]
        end
    end
    fill!(weights, 1 / Np)
    return weights, particles
end
