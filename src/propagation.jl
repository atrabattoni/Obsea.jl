function tlim(incidence, mode, depth, celerity)
    (2 * mode - 1) * depth / cos(incidence) / celerity
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


function tdoa(r, mode, depth, celerity, ic, ib)
    toa(r, mode + 1, depth, celerity, ic, ib) -
    toa(r, mode, depth, celerity, ic, ib)
end

# struct TdoaLut
#     v::Dict{Int64,Array{Float64,1}}
#     τ::Array{Float64,2}
#     function TdoaTable(propagation)
#         v = Dict(mode => window(6*σ[mode], σ[mode]) for mode = 1:length(n))
#         τ = [
#             tdoa(r, mode, depth, celerity, ic, ib)
#             for r in rrange, mode = 1:nmode
#         ]
#         new(nmode, depth, celerity, ic, ib, lam, v)
#     end
# end
