import Obsea: State, Grid
import Plots.plot
using StructArrays

function Plots.plot(cloud::StructVector{State}, grid::Grid)
    rmax = last(grid.r)
    mask = (cloud.m .!= 0)
    scatter(
        cloud[mask].x / 1000.0,
        cloud[mask].y / 1000.0;
        xlim = (-rmax / 1000.0, rmax / 1000.0),
        ylim = (-rmax / 1000.0, rmax / 1000.0),
        alpha = 0.5,
        aspect_ratio = :equal,
    )
end

function Plots.plot!(cloud::StructVector{State}, grid::Grid)
    rmax = last(grid.r)
    mask = (cloud.m .!= 0)
    scatter!(
        cloud[mask].x / 1000.0,
        cloud[mask].y / 1000.0;
        xlim = (-rmax / 1000.0, rmax / 1000.0),
        ylim = (-rmax / 1000.0, rmax / 1000.0),
        alpha = 0.5,
        aspect_ratio = :equal,
    )
end

function Plots.plot!(cloud::StructVector{State}, grid::Grid, c)
    rmax = last(grid.r)
    mask = (cloud.m .!= 0)
    scatter!(
        cloud[mask].x / 1000.0,
        cloud[mask].y / 1000.0;
        zcolor = c,
        xlim = (-rmax / 1000.0, rmax / 1000.0),
        ylim = (-rmax / 1000.0, rmax / 1000.0),
        alpha = 0.5,
        aspect_ratio = :equal,
    )
end
