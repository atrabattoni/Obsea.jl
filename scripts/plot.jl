import Obsea: State, Grid
using Plots
using StructArrays

function Plots.scatter!(x, y, grid::Grid; kw...)
    rmax = last(grid.r) / 1000
    x = x / 1000
    y = y / 1000
    scatter!(
        x,
        y;
        xlim = (-rmax, rmax),
        ylim = (-rmax, rmax),
        aspect_ratio = :equal,
        markerstrokealpha = 0,
        kw...,
    )
end

function Plots.scatter!(cloud::StructVector{State}, grid::Grid; kw...)
    rmax = last(grid.r) / 1000
    x = cloud.x / 1000
    y = cloud.y / 1000
    mask = (cloud.m .!= 0)
    scatter!(
        x,
        y;
        xlim = (-rmax, rmax),
        ylim = (-rmax, rmax),
        aspect_ratio = :equal,
        alpha = 0.5,
        kw...,
    )
end

function Plots.plot!(
    particles::StructArray{State,2},
    grid::Grid;
    N = 100,
    kw...,
)
    rmax = last(grid.r) / 1000
    Np = size(particles, 2)
    idx = rand(1:Np, N)
    x = replace!(particles.x[:, idx], 0.0 => NaN) ./ 1000
    y = replace!(particles.y[:, idx], 0.0 => NaN) ./ 1000
    plot!(
        x,
        y;
        xlim = (-rmax, rmax),
        ylim = (-rmax, rmax),
        aspect_ratio = :equal,
        legend = false,
        kw...,
    )
end


function Plots.heatmap!(z::Array{Float64,2}, xgrid::StepRangeLen; kw...)
    rmax = maximum(abs.(collect(xgrid))) / 1000
    x = y = xgrid / 1000
    z = log.(abs.(z))
    heatmap!(
        x,
        y,
        z;
        xlim = (-rmax, rmax),
        ylim = (-rmax, rmax),
        aspect_ratio = :equal,
        clim = (-20, 20),
        kw...,
    )
end
