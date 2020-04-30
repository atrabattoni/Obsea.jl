function update!(weights, cloud, ℓ, models, grid)
    pd = [model.pd for model in models]
    ℓr, ℓa = make(ℓ, grid)
    r = sqrt.(cloud.x.^2 .+ cloud.y.^2)
    a = mod.(atan.(cloud.x, cloud.y), 2π)
    weights .*= (1.0 .- pd) .+ pd .* ℓr.(r, cloud.m) .* ℓa.(cloud.f, a, cloud.m)
    weights ./= sum(weights)
end

function make(ℓ, grid)
    ℓr = interpolate!(ℓ.r, (BSpline(Linear()), NoInterp()))
    ℓa = interpolate!(ℓ.a, (BSpline(Linear()), BSpline(Linear()), NoInterp()))
    ℓr = extrapolate(ℓr, 1.0)
    ℓa = extrapolate(ℓa, 1.0)
    ℓr = scale(ℓr, grid2range(grid.r), 1:grid.Nm)
    ℓa = scale(ℓa, grid2range(grid.f), grid2range(grid.a), 1:grid.Nm)
    ℓr, ℓa
end
