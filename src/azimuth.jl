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


struct Azimuth
    mrl::Float64
    n::Int64
    arange::AbstractRange
end


function precomp(z, model::Azimuth)
    @unpack mrl, n, arange = model
    y0 = [wrapcauchy(z[i], a, mrl) for i = 1:length(z), a in arange]
    y1 = rollprod(y0, n)
    cat(y0, y1, dims = 3)
end
