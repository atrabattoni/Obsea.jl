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

function convsame7(u, v)
    out = zeros(length(u))
    @inbounds for j = 4:length(out)-3
        out[j] += u[j-3] * v[1]
        out[j] += u[j-2] * v[2]
        out[j] += u[j-1] * v[3]
        out[j] += u[j] * v[4]
        out[j] += u[j+1] * v[5]
        out[j] += u[j+2] * v[6]
        out[j] += u[j+3] * v[7]
    end
    out
end


function convsame8(u, v)
    pad = length(v) ÷ 2
    out = zeros(length(u))
    @avx for i = 1:length(v)
        for j = 1+pad:length(out)-pad
            out[j] += u[j+i-pad-1] * v[i]
        end
    end
    out
end

function convsame9(u, v)
    pad = length(v) ÷ 2
    out = Vector{Float64}(undef, length(u))
    @avx for j = 1+pad:length(out)-pad
        s = 0.0
        for i = 1:length(v)
            s += u[j+i-pad-1] * v[i]
        end
        out[j] = s
    end
    out
end


@assert convsame2(u, v) ≈ convsame1(u, v)
@assert convsame3(u, v) ≈ convsame1(u, v)
@assert convsame4(u, v) ≈ convsame1(u, v)
@assert convsame5(u, v) ≈ convsame1(u, v)
@assert convsame6(u, v)[1+3:end-3] ≈ convsame1(u, v)[1+3:end-3]
@assert convsame7(u, v)[1+3:end-3] ≈ convsame1(u, v)[1+3:end-3]
@assert convsame8(u, v)[1+3:end-3] ≈ convsame1(u, v)[1+3:end-3]
@assert convsame9(u, v)[1+3:end-3] ≈ convsame1(u, v)[1+3:end-3]

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
println("avx no padding")
@btime convsame6(u, v);
println("inbounds unrolling")
@btime convsame7(u, v);
println("avx other order looping")
@btime convsame8(u, v);
println("avx accumulator")
@btime convsame9(u, v);
println()

@profiler for i = 1:10000
    convsame9(u, v)
end
