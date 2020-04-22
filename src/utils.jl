function window(n::Int64, σ::Float64)
    v = gaussian(n, σ)
    v /= sum(v)
end


function limit(xrange, xmin, xmax)
    imin = xmin ÷ xrange.step + 1
    imax = xmax ÷ xrange.step
    xrange[imin:imax]
end


function symbolize(d)
    Dict(Symbol(k) => v for (k, v) in d)
end


function convsame(u, v)
    @assert length(u) > length(v)
    @assert isodd(length(v))
    padding = length(v) ÷ 2
    conv(u, v)[1+padding:end-padding]
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


function argsample(cdf)
    searchsortedfirst(cdf, rand(), lt = <=)
end


function argsample(cdf, N)
    out = Array{Int64,1}(undef, N)
    j = 1
    s = cdf[j]
    u = rand() / N
    for i = 1:N
        while s < u
            j += 1
            s += cdf[j]
        end
        out[i] = j
        u += 1 / N
    end
    out
end
