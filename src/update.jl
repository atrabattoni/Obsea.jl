function Scan(y_r, y_fam, grid)
    itp_r = interpolate(y_r, BSpline(Cubic(Line(OnGrid()))))
    itp_fam = interpolate(
        y_fam,
        (
            BSpline(Cubic(Line(OnGrid()))),
            BSpline(Cubic(Periodic(OnGrid()))),
            NoInterp(),
        ),
    )
    itp_r = extrapolate(itp_r, 1.0)
    itp_fam = extrapolate(itp_fam, (1.0, Periodic(), Throw()))
    itp_r = scale(itp_r, grid.range_r)
    itp_fam = scale(itp_fam, grid.range_f, grid.range_a, grid.range_m)

    cdf_r = cumsum(y_r)
    cdf_fam = cumsum(vec(y_fam))
    cdf_r /= cdf_r[end]
    cdf_fam /= cdf_fam[end]

    Scan(cdf_r, cdf_fam, itp_r, itp_fam)
end


function update!(weights, cloud, scan, params)
    for (i, particle) in enumerate(cloud)
        weights[i] *= exp(logl(scan, particle[end], params))
    end
    weights /= sum(weights)
end


function logl(scan, state, params)
    @unpack pd = params
    if isempty(state)
        return 0.0
    else
        r = sqrt(state.x^2 + state.y^2)
        f = state.frequency
        a = atan(state.x, state.y)
        if state.model == 1
            return log(scan.itp_r(r) * scan.itp_fam(f, a, 1))
        elseif state.model == 2
            return log((1.0 - pd) + pd * scan.itp_r(r) * scan.itp_fam(f, a, 2))
        end
    end
end
