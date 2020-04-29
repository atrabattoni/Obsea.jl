using BenchmarkTools
using Profile

function argsample1(ℓ, N)
    scale = sum(ℓ)
    step = scale / N
    u = rand() * step
    out = Vector{Int}(undef, N)
    j = 1
    s = ℓ[j]
    for i = 1:N
        while s < u
            j += 1
            s += ℓ[j]
        end
        out[i] = j
        u += step
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


function argsample2(ℓ, N)
    F = cumsum(ℓ)
    scale = last(F)
    out = Vector{Int}(undef, N)
    j = 1
    u = rand() * scale / N
    s = ℓ[j]
    for i = 1:N
        while s < u
            j += 1
            s = F[j]
        end
        out[i] = j
        u += scale / N
    end
    out
end

function argsample3(F, N)
    scale = last(F)
    out = Vector{Int}(undef, N)
    j = 1
    u = rand() * scale / N
    s = ℓ[j]
    for i = 1:N
        while s < u
            j += 1
            s = F[j]
        end
        out[i] = j
        u += scale / N
    end
    out
end

N = 100000
ℓ = ones(N)
F = cumsum(ℓ)
scale = float(N)
@assert(argsample1(ℓ, N) == collect(1:N))
@assert(argsample(ℓ, N) == collect(1:N))
@assert(argsample2(ℓ, N) == collect(1:N))
@assert(argsample3(F, N) == collect(1:N))

println()
println("reference")
@btime argsample1($ℓ, $N);
println("reference with scale")
@btime argsample($ℓ, $N);
println("cumsum")
@btime argsample2($ℓ, $N);
println("precomputed cumsum")
@btime argsample3($F, $N);
println()

@profiler for i = 1:1000
    argsample(ℓ, N)
end

@assert argsample(ℓ, N, scale = 2 * sum(ℓ))[N÷2+1:end] == fill(0, N ÷ 2)
