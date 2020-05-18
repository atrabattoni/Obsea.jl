function update!(weights, cloud, ℓt, models, grid)
    @unpack pd = models
    ℓr, ℓa = make(ℓt, grid)
    mask = .!isdead.(cloud)
    @views weights = weights[mask]
    @views cloud = cloud[mask]
    r = sqrt.(cloud.x .^ 2 .+ cloud.y .^ 2)
    f = cloud.f
    a = mod.(atan.(cloud.x, cloud.y), 2π)
    m = cloud.m
    weights .*= (1.0 .- pd[m]) .+ pd[m] .* ℓr.(r, m) .* ℓa.(f, a, m)
    weights ./= sum(weights)
end

function make(ℓt, grid)
    ℓr = interpolate!(ℓt.r, (BSpline(Linear()), NoInterp()))
    ℓa = interpolate!(ℓt.a, (BSpline(Linear()), BSpline(Linear()), NoInterp()))
    ℓr = extrapolate(ℓr, 1.0)
    ℓa = extrapolate(ℓa, 1.0)
    ℓr = scale(ℓr, grid2range(grid.r), 1:grid.Nm)
    ℓa = scale(ℓa, grid2range(grid.f), grid2range(grid.a), 1:grid.Nm)
    ℓr, ℓa
end
