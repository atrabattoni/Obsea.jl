function window(sigma)
    if sigma == 0.0
        return [1.0]
    else
        n = 2 * round(Int, 2 * sigma) + 1
        v = gaussian(n, sigma / n)
        v /= sum(v)
        return v
    end
end


function limit(grid, min, max)
    @assert first(grid) <= min <= max <= last(grid)
    i = 1
    for x in grid
        if x >= min
            break
        else
            i += 1
        end
    end
    j = length(grid)
    for x in reverse(grid)
        if x <= max
            break
        else
            j -= 1
        end
    end
    grid[i:j]
end


function grid2range(grid)
    range(first(grid), last(grid), length = length(grid))
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
    @inbounds for j = Nc:Nu-Nc
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

function argsample(ℓ, N; scale = sum(ℓ))
    step = scale / N
    u = rand() * step
    out = Vector{Int}(undef, N)
    j = 1
    s = ℓ[j]
    for i = 1:N
        while s < u
            j += 1
            if j <= length(ℓ)
                s += ℓ[j]
            else
                s = Inf
                j = 0
            end
        end
        out[i] = j
        u += step
    end
    out
end

function ra2xy(r, a)
    r * sin(a), r * cos(a)
end

function xy2ra(x, y)
    sqrt(x^2 + y^2), atan(x, y) % 2π
end
