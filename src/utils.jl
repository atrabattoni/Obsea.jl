function window(n, sigma)
    v = gaussian(n, sigma)
    v /= sum(v)
end


function limit(xrange, xmin, xmax)
    @assert first(xrange) <= xmin <= xmax <= last(xrange)
    imin = 1
    for x in xrange
        if x >= xmin
            break
        else
            imin += 1
        end
    end
    imax = length(xrange)
    for x in reverse(xrange)
        if x <= xmax
            break
        else
            imax -= 1
        end
    end
    xrange[imin:imax]
end


function symbolize(d)
    Dict(Symbol(k) => v for (k, v) in d)
end


function convsame(u, v)
    Nu = length(u)
    Nv = length(v)
    @assert Nu > Nv
    @assert isodd(Nv)
    Nc = (Nv ÷ 2) + 1
    out = similar(u)
    @inbounds for j = 1:Nc-1
        s = 0.0
        for i = Nc-(j-1):Nv
            s += u[j+i-Nc] * v[i]
        end
        out[j] = s
    end
    @avx for j = Nc:Nu-Nc
        s = 0.0
        for i = 1:Nv
            s += u[j+i-Nc] * v[i]
        end
        out[j] = s
    end
    @inbounds for j = Nu-Nc+1:Nu
        s = 0.0
        for i = 1:Nc+(Nu-j)
            s += u[j+i-Nc] * v[i]
        end
        out[j] = s
    end
    return out
end


function wrapcauchy(z, a, mrl)
    @assert 0 < mrl < 1.0
    (1.0 - mrl^2) / (1.0 - 2.0 * mrl * cos(z - a) + mrl^2)
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


function argsample(cdf; scale = 1.0)
    @assert 0.0 <= first(cdf) <= last(cdf) <= scale + 10 * eps(scale)
    rng = rand() * scale
    if rng > last(cdf)
        return 0
    else
        return searchsortedfirst(cdf, rng, lt = <=)
    end
end


function argsample(cdf, N; scale = 1.0)
    @assert 0.0 <= first(cdf) <= last(cdf) <= scale + 10 * eps(scale)
    out = Array{Int64,1}(undef, N)
    j = 1
    s = cdf[j]
    u = rand() * scale / N
    for i = 1:N
        while (j != 0) & (s < u)
            j += 1
            if !(j > length(cdf))
                s += cdf[j]
            else
                j = 0
            end
        end
        out[i] = j
        u += scale / N
    end
    out
end

function ra2xy(r, a)
    r * sin(a), r * cos(a)
end

function xy2ra(x, y)
    sqrt(x^2 + y^2), atan(x, y) % 2π
end
