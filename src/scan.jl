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
