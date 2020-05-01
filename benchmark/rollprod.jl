using BenchmarkTools
using LoopVectorization

function rollprod1(out, u, Nv)
    @assert isodd(Nv)
    for i = 1:length(u)
        imin = max(1, i - (Nv ÷ 2))
        imax = min(size(u, 1), i + (Nv ÷ 2))
        @views out[i] = prod(u[imin:imax])
    end
    return out
end

function rollprod2(out, u, Nv)
    Nu = length(u)
    @assert isodd(Nv)
    Np = Nv ÷ 2
    for j = 1:Np
        @views out[j] = prod(u[1:j+Np])
    end
    for j = Np+1:Nu-Np-1
        @views out[j] = prod(u[j-Np:j+Np])
    end
    for j = Nu-Np:Nu
        @views out[j] = prod(u[j-Np:Nu])
    end
    return out
end

function rollprod3(out, u, Nv)
    Nu = length(u)
    @assert isodd(Nv)
    Np = Nv ÷ 2
    @inbounds for j = 1:Np
        p = 1.0
        for i = 1:j+Np
            p *= u[i]
        end
        out[j] = p
    end
    @avx for j = Np+1:Nu-Np-1
        p = 1.0
        for i = -Np:Np
            p *= u[j+i]
        end
        out[j] = p
    end
    @inbounds for j = Nu-Np:Nu
        p = 1.0
        for i = j-Np:Nu
            p *= u[i]
        end
        out[j] = p
    end
    return out
end

function rollprod4(out, u, Nv)
    Nu = length(u)
    @assert isodd(Nv)
    Nc = (Nv ÷ 2) + 1
    out .= log.(u)
    cum = cumsum(out)
    for j = 1:Nc
        out[j] = cum[j + Nc - 1]
    end
    for j = Nc+1:Nu-Nc
        out[j] = cum[j + Nc - 1] - cum[j - Nc]
    end
    for j = Nu-Nc+1:Nu
        out[j] = cum[Nu] - cum[j - Nc]
    end
    out .= exp.(out)
    return out
end

u = rand(10_000_001)
out = similar(u)
Nv = 5

@assert all(rollprod2(out, u, Nv) .≈ rollprod1(out, u, Nv))
@assert all(rollprod3(out, u, Nv) .≈ rollprod1(out, u, Nv))
@assert all(rollprod4(out, u, Nv) .≈ rollprod1(out, u, Nv))


println()
println("reference")
@btime rollprod1($out, $u, $Nv)
println("splited loop")
@btime rollprod2($out, $u, $Nv)
println("fully looped aux petits onions")
@btime rollprod3($out, $u, $Nv)
println("log cumsum")
@btime rollprod4($out, $u, $Nv)
println()
