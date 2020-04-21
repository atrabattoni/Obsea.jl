function tlim(incidence, mode, depth, celerity)
    (2 * mode - 1) * depth / cos(incidence) / celerity
end


function toa(r, mode, depth, celerity)
    sqrt(((2 * mode - 1) * depth)^2 + r^2) / celerity
end

function toa(r, mode, depth, celerity, ic, ib)
    t = toa(r, mode, depth, celerity)
    tc = tlim(ic, mode, depth, celerity)
    tb = tlim(ib, mode, depth, celerity)
    if tc <= t <= tb
        return t
    else
        return NaN
    end
end


function tdoa(r, mode, depth, celerity, ic, ib)
    toa(r, mode + 1, depth, celerity, ic, ib) -
    toa(r, mode, depth, celerity, ic, ib)
end


function propagation(rrange, nmode, depth, celerity, ic, ib)
    [tdoa(r, mode, depth, celerity, ic, ib) for r in rrange, mode = 1:nmode]
end


function window(n::Int64, σ::Float64)
    v = gaussian(n, σ)
    v /= sum(v)
end

function window(n::Array{Int64,1}, σ::Array{Float64,1})
    @assert length(n) == length(σ)
    Dict(mode => window(n[mode], σ[mode]) for mode = 1:length(n))
end


function convsame(u, v)
    @assert length(u) > length(v)
    @assert isodd(length(v))
    padding = length(v) ÷ 2
    conv(u, v)[1+padding:end-padding]
end


struct Range
    nmode::Int64
    depth::Float64
    celerity::Float64
    ic::Float64
    ib::Float64
    lam::Float64
    n::Array{Float64,1}
    σ::Array{Int64,1}
    rrange::AbstractRange
    τrange::AbstractRange
    v::Dict{Int64,Array{Float64,1}}
    τ::Array{Float64,2}
    function Range(nmode, depth, celerity, ic, ib, lam, n, σ, rrange, τrange)
        @assert length(n) == length(σ) == nmode
        v = window(n, σ)
        τ = propagation(rrange, nmode, depth, celerity, ic, ib)
        new(nmode, depth, celerity, ic, ib, lam, n, σ, rrange, τrange, v, τ)
    end
end


function precomp(z::Array{Float64,1}, model::Range)
    @unpack lam, nmode, v, τrange, τ = model
    @assert size(τ, 2) == length(v) == nmode
    u = exp.(-lam ./ 2.0) .* besseli.(0, sqrt.(lam .* z))
    y = ones(size(τ, 1))
    for mode = 1:nmode
        A = convsame(u, v[mode])
        itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
        itp = extrapolate(itp, 1.0)
        itp = scale(itp, τrange)
        y .*= itp(τ[:, mode])
    end
    y
end
