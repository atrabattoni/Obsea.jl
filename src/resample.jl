function argsample(weights, N)
    out = Array{Int64,1}(undef, N)
    j = 1
    s = weights[j]
    u = rand() / N
    for i = 1:N
        while s < u
            j += 1
            s += weights[j]
        end
        out[i] = j
        u += 1 / N
    end
    out
end


function resample(weights, cloud)
    @assert length(weights) == length(cloud)
    idx = argsample(weights, length(cloud))
    newcloud = [deepcopy(cloud[i]) for i in idx]
    newweights = fill(1 / length(cloud), length(cloud))
    (newweights, newcloud)
end
