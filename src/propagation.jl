struct TDOALUT
    Nmode::Int
    v::Dict{Int,Vector{Float64}}
    τ::Array{Float64,2}
    function TDOALUT(propa, grid)
        @unpack Nmode, depth, celerity, ic, ib, sigma = propa
        v = Dict(mode => window(sigma[mode]) for mode = 1:Nmode)
        τ = [
            tdoa(r, mode, depth, celerity, ic, ib)
            for r in grid.r, mode = 1:Nmode
        ]
        new(Nmode, v, τ)
    end
end

function tdoa(r, mode, depth, celerity, ic, ib)
    toa(r, mode + 1, depth, celerity, ic, ib) -
    toa(r, mode, depth, celerity, ic, ib)
end

function toa(r, mode, depth, celerity, ic, ib)
    t = sqrt(((2 * mode - 1) * depth)^2 + r^2) / celerity
    tc = tlim(ic, mode, depth, celerity)
    tb = tlim(ib, mode, depth, celerity)
    if tc <= t <= tb
        return t
    else
        return NaN
    end
end

function tlim(incidence, mode, depth, celerity)
    (2 * mode - 1) * depth / cos(incidence) / celerity
end
