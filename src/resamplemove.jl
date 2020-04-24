function resample(weights, cloud)
    @assert length(weights) == length(cloud)
    cdf = cumsum(weights)
    idx = argsample(cdf, length(cloud))
    newcloud = [deepcopy(cloud[i]) for i in idx]
    newweights = fill(1 / length(cloud), length(cloud))
    newweights, newcloud
end
