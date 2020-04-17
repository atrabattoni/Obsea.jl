function argsample(w, N)
    out = Array{Int64,1}(undef, N)
    j = 1
    s = w[j]
    u = rand() / N
    for i=1:N
        while s<u
            j += 1
            s += w[j]
        end
        out[i] = j
        u += 1 / N
    end
    out
end


function resample(cloud)
    w = [particle.metadata.weight for particle ∈ cloud]
    idx = argsample(w, length(cloud))
    out = Cloud()
    for i ∈ idx
        particle = deepcopy(cloud[i])
        particle.metadata.weight = 1 / length(cloud)
        push!(out, particle)
    end
    out
end
