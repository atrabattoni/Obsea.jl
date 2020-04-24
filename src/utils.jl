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
    @assert length(u) > length(v)
    @assert isodd(length(v))
    pad = length(v) ÷ 2
    out = zeros(length(u))
    u = [zeros(pad); u; zeros(pad)]
    @avx for j = 1:length(out)
        for i = 1:length(v)
            out[j] += u[j+i-1] * v[i]
        end
    end
    out
end


function wrapcauchy(z, a, mrl)
    @assert 0 < mrl < 1.0
    (1.0 - mrl^2) / (1.0 - 2.0 * mrl * cos(z - a) + mrl^2)
end


function rollprod(x, n)
    @assert isodd(n)
    out = similar(x)
    for j = 1:size(x, 2)
        for i = 1:size(x, 1)
            imin = max(1, i - (n ÷ 2))
            imax = min(size(x, 1), i + (n ÷ 2))
            out[i, j] = prod(x[imin:imax, j])
        end
    end
    out
end


function argsample(cdf; scale=1.0)
    @assert 0.0 <= first(cdf) <= last(cdf) <= scale + 10 * eps(scale)
    rng = rand() * scale
    if rng > last(cdf)
        return 0
    else
        return searchsortedfirst(cdf, rng, lt = <=)
    end
end


function argsample(cdf, N; scale=1.0)
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
