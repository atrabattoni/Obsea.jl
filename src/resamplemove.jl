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
            k = pop!(stack)
            @views particles[:, j] .= particles[:, k]
        end
    end
    fill!(weights, 1 / Np)
    return weights, particles
end
