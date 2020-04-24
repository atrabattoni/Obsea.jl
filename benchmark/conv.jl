using BenchmarkTools
using DSP
using LoopVectorization

u = rand(513)
v = gaussian(7, 1 / 6)

function convsame1(u, v)
    pad = length(v) ÷ 2
    conv(u, v)[1+pad:end-pad]
end

function convsame2(u, v)
    pad = length(v) ÷ 2
    out = zeros(length(u))
    u = [zeros(pad); u; zeros(pad)]
    for j = 1:length(out)
        for i = 1:length(v)
            out[j] += u[j+i-1] * v[i]
        end
    end
    out
end

function convsame3(u, v)
    pad = length(v) ÷ 2
    out = zeros(length(u))
    u = [zeros(pad); u; zeros(pad)]
    @inbounds @simd for j = 1:length(out)
        for i = 1:length(v)
            out[j] += u[j+i-1] * v[i]
        end
    end
    out
end

function convsame4(u, v)
    pad = length(v) ÷ 2
    out = zeros(length(u))
    u = [zeros(pad); u; zeros(pad)]
    @inbounds for j = 1:length(out)
        for i = 1:length(v)
            out[j] += u[j+i-1] * v[i]
        end
    end
    out
end

function convsame5(u, v)
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

function convsame6(u, v)
    pad = length(v) ÷ 2
    out = zeros(length(u))
    @avx for j = 1+pad:length(out)-pad
        for i = 1:length(v)
            out[j] += u[j+i-pad-1] * v[i]
        end
    end
    out
end

@assert convsame2(u, v) ≈ convsame1(u, v)
@assert convsame3(u, v) ≈ convsame1(u, v)
@assert convsame4(u, v) ≈ convsame1(u, v)
@assert convsame5(u, v) ≈ convsame1(u, v)
@assert convsame6(u, v)[1+3:end-3] ≈ convsame1(u, v)[1+3:end-3]

println()
println("reference")
@btime convsame1(u, v);
println("simple implementation")
@btime convsame2(u, v);
println("inbounds")
@btime convsame3(u, v);
println("inbounds simd")
@btime convsame4(u, v);
println("avx")
@btime convsame5(u, v);
println("no padding")
@btime convsame6(u, v);
println()

@profiler for i = 1:10000
    convsame6(u, v)
end
