function resample!(weights, particles, ℓ, models, grid)
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
            @views mcmc!(particles[:, j], particles[:, k], ℓ, models, grid)
        end
    end
    fill!(weights, 1 / Np)
    return weights, particles
end

function mcmc!(particle, prevparticle, ℓ, models, grid)
    particle .= prevparticle
end
