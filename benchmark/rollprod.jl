using BenchmarkTools

function rollprod1(u, Nv)
    @assert isodd(Nv)
    out = similar(u)
    for i = 1:length(u)
        imin = max(1, i - (Nv ÷ 2))
        imax = min(size(u, 1), i + (Nv ÷ 2))
        out[i] = prod(u[imin:imax])
    end
    out
end

function rollprod2(u, Nv)
    Nu = length(u)
    @assert isodd(Nv)
    Nc = (Nv ÷ 2) + 1
    out = similar(u)
    cum = cumsum(log.(u))
    for j = 1:Nc
        out[j] = cum[j + Nc - 1]
    end
    for j = Nc+1:Nu-Nc
        out[j] = cum[j + Nc - 1] - cum[j - Nc]
    end
    for j = Nu-Nc+1:Nu
        out[j] = cum[Nu] - cum[j - Nc]
    end
    exp.(out)
end

function rollprod(u, Nv)
    Nu = length(u)
    @assert isodd(Nv)
    Np = Nv ÷ 2
    out = similar(u)
    for j = 1:Np
        out[j] = prod(u[1:j+Np])
    end
    for j = Np+1:Nu-Np-1
        out[j] = prod(u[j-Np:j+Np])
    end
    for j = Nu-Np:Nu
        out[j] = prod(u[j-Np:Nu])
    end
    out
end


u = rand(513)
Nv = 7

@assert all(rollprod2(u, Nv) .≈ rollprod1(u, Nv))
@assert rollprod(u, Nv) == rollprod1(u, Nv)

println()
println("reference")
@btime rollprod1($u, $Nv)
println("splited loop")
@btime rollprod($u, $Nv)
println("cumsum")
@btime rollprod2($u, $Nv)
println()
