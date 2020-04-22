function resample!(weights, cloud)
    @assert length(weights) == length(cloud)
    idx = argsample(weights, length(cloud))
    newcloud = [deepcopy(cloud[i]) for i in idx]
    newweights = fill(1 / length(cloud), length(cloud))
    fill!(weights, newweights)
    fill!(cloud, newcloud)
end
