struct State
    model::Int64
    frequency::Float64
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end


struct EmptyState end

Trajectory = Array{Union{State,EmptyState},1}


mutable struct Metadata
    weight::Float64
end


struct Particle
    trajectory::Trajectory
    metadata::Metadata
end


struct Parameters
    T::Float64
    q::Float64
    ps::Float64
    pb::Float64
end


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
        itp_r = scale(itp_r, grid.range_r)
        itp_fam = scale(itp_fam, grid.range_f, grid.range_a, grid.range_m)
        # itp_r = extrapolate(itp_r, 1.0)
        # itp_fam = extrapolate(itp_fam, (1.0, Periodic(), Throw()))
        new(itp_r, itp_fam)
    end
end
