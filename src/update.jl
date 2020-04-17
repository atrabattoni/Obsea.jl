struct Grid
    range_r::AbstractRange
    range_f::AbstractRange
    range_a::AbstractRange
    range_m::UnitRange{Int64}
end


struct Scan
    itp_r::ScaledInterpolation
    itp_fam::ScaledInterpolation
    function Scan(y_r::Array{Float64,1}, y_fam::Array{Float64,3}, grid::Grid)
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
        new(itp_r, itp_fam)
    end
end


function update!(cloud::Cloud, scan::Scan, params::Parameters)
    normalization = 0.0
    for particle ∈ cloud
        likelihood = exp(logl(scan, particle.trajectory[end], params))
        particle.metadata.weight *= likelihood
        normalization += likelihood
    end
    for particle ∈ cloud
        particle.metadata.weight /= normalization
    end
end


function logl(scan::Scan, state::ShipState, params::Parameters)
    pd = params.pd
    r = sqrt(state.x^2 + state.y^2)
    f = state.frequency
    a = atan(state.x, state.y)
    log(scan.itp_r(r) * scan.itp_fam(f, a, 1))
end

function logl(scan::Scan, state::WhaleState, params::Parameters)
    pd = params.pd
    r = sqrt(state.x^2 + state.y^2)
    f = state.frequency
    a = atan(state.x, state.y)
    log((1.0 - pd) + pd * scan.itp_r(r) * scan.itp_fam(f, a, 2))
end

function logl(scan::Scan, state::EmptyState)
    0.0
end
