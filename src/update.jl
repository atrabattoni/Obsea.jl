function update!(weights, cloud, ℓt, models, grid)
    ℓ = make(ℓt, models, grid)
    mask = .!isdead.(cloud)
    @views cloud = cloud[mask]
    r = sqrt.(cloud.x .^ 2 .+ cloud.y .^ 2)
    f = cloud.f
    a = mod.(atan.(cloud.x, cloud.y), 2π)
    m = cloud.m
    weights[mask] .*= ℓ.(r, f, a, m)
    weights ./= sum(weights)
end

function make(ℓt, models, grid)
    @unpack pd = models
    ℓr = interpolate!(ℓt.r, (BSpline(Linear()), NoInterp()))
    ℓa = interpolate!(ℓt.a, (BSpline(Linear()), BSpline(Linear()), NoInterp()))
    ℓr = extrapolate(ℓr, 1.0)
    ℓa = extrapolate(ℓa, 1.0)
    ℓr = scale(ℓr, grid2range(grid.r), 1:grid.Nm)
    ℓa = scale(ℓa, grid2range(grid.f), grid2range(grid.a), 1:grid.Nm)
    return ℓ(r, f, a, m) = (1.0 .- pd[m]) .+ pd[m] .* ℓr.(r, m) .* ℓa.(f, a, m)
end
