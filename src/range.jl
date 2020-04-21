function tlim(incidence, mode, celerity, depth)
    (2 * mode - 1) * depth / cos(incidence) / celerity
end


function toa(r, mode, celerity, depth, ic, ib)
    t = sqrt((mode * depth)^2 + r^2) / celerity
    tc = tlim(ic, mode, celerity, depth)
    tb = tlim(ib, mode, celerity, depth)
    if tc < t < tb
        return t
    else
        return NaN
    end
end


function tdoa(r, mode, celerity, depth, ic, ib)
    toa(r, mode + 1, celerity, depth, ic, ib) -
    toa(r, mode, celerity, depth, ic, ib)
end


function propagation(rrange, nmode, celerity, depth, ic, ib)
    [tdoa(r, mode, celerity, depth, ic, ib) for r in rrange, mode = 1:nmode]
end


function window(n::Int64, σ::Float64)
    v = gaussian(n[mode], σ[mode])
    v /= sum(v)
end

function window(n::Array{Int64,1}, σ::Array{Float64,1})
    @assert length(n) = length(σ)
    Dict(mode => window(n[mode], σ[mode]) for mode = 1:length(n))
end


function convsame(u, v)
    @assert length(u) > length(v)
    @assert isodd(length(v))
    padding = length(v) ÷ 2
    conv(u, v)[1+padding:end-padding]
end


function precomp(z, lam, nmode, v, τrange, τ)
    u = exp.(-lam ./ 2.0) .* besseli.(0, sqrt.(lam .* z))
    y = ones(size(τ, 1))
    for mode = 1:nmode
        A = convsame(u, v[mode])
        itp = interpolate!(A, BSpline(Cubic(Line(OnGrid()))))
        itp = extrapolate(itp, 1.0)
        itp = scale(itp, τrange)
        y .*= itp(τ[:, mode])
    end
end


struct Range
    nmode::Int64
    celerity::Float64
    depth::Float64
    ic::Float64
    ib::Float64
    lam::Float64
    n::Array{Float64,1}
    σ::Array{Int64,1}
    rrange::AbstractRange
    τrange::AbstractRange
    v::Dict{Int64,Array{Float64,1}}
    τ::Array{::Float64,2}
    function Range(nmode, celerity, depth, ic, ib, lam, n, σ, rrange, τrange)
        v = windows(n, σ)
        τ = propagation(rrange, nmode, celerity, depth, ic, ib)
        new(nmode, celerity, depth, ic, ib, lam, n, σ, rrange, τrange, v, τ)
    end
end


function precomp(z, model::Range)
    @unpack z, lam, nmode, n, σ, τrange, τ = model
    precomp(z, lam, nmode, n, σ, τrange, τ)
end
