function resample(weights, cloud)
    idx = argsample(weights, length(cloud); scale = 1.0)
    newcloud = [deepcopy(cloud[i]) for i in idx]
    newweights = fill(1 / length(cloud), length(cloud))
    newweights, newcloud
end
