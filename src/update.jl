function update!(weights, cloud, ℓ, models, grid)
    for (i, particle) in enumerate(cloud)
        weights[i] *= interp(ℓ, last(particle), models, grid)
    end
    weights ./= sum(weights)
end


function interp(ℓ, state, models, grid)
    if isempty(state)
        return 1.0
    else
        @unpack m, f, x, y = state
        @unpack pd = models[m]
        @unpack Nm = grid
        r, a = xy2ra(x, y)

        ℓr = interpolate!(ℓ.r, (BSpline(Linear()), NoInterp()))
        ℓr = extrapolate(ℓr, 1.0)
        ℓr = scale(ℓr, grid2range(grid.r), 1:Nm)

        ℓa = interpolate!(
            ℓ.a,
            (BSpline(Linear()), BSpline(Linear()), NoInterp()),
        )
        ℓa = extrapolate(ℓa, 1.0)
        ℓa = scale(ℓa, grid2range(grid.f), grid2range(grid.a), 1:Nm)

        return (1.0 - pd) + pd * ℓr(r, m) * ℓa(f, a, m)
    end
end
