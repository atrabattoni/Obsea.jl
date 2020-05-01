function resample!(weights, particles)
    idx = argsample(weights, length(weights))
    Np = length(weights)
    keep = zeros(Bool, Np)
    stack = Int[]
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
