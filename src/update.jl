function update!(weights, cloud, itp, models)
    for (i, particle) in enumerate(cloud)
        weights[i] *= likelihoodratio(itp, last(particle), models)
    end
    weights ./= sum(weights)
end


function likelihoodratio(itp, state, models)
    if isempty(state)
        return 1.0
    else
        @unpack m, f, x, y = state
        @unpack pd = models[m]
        r, a = xy2ra(x, y)
        return (1.0 - pd) + pd * itp(r, f, a, m)
    end
end
