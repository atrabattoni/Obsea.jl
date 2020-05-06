using BenchmarkTools
using StaticArrays
using LoopVectorization

function wrapcauchy(za, a, mrl)
    Nf, Nt = size(za)
    Na = length(a)
    Nm = length(mrl)

    za = reshape(za, (Nf, 1, Nt, 1))
    a = reshape(a, (1, Na, 1, 1))
    mrl = reshape(mrl, (1, 1, 1, Nm))

    return (1.0 .- mrl .^ 2) ./ (1.0 .- 2.0 .* mrl .* cos.(za .- a) .+ mrl .^ 2)
end

function wrapcauchy!(out, za, a, mrl)
    Nf, Nt = size(za)
    Na = length(a)
    Nm = length(mrl)

    za = reshape(za, (Nf, 1, Nt, 1))
    a = reshape(a, (1, Na, 1, 1))
    mrl = reshape(mrl, (1, 1, 1, Nm))

    out .= (1.0 .- mrl .^ 2) ./ (1.0 .- 2.0 .* mrl .* cos.(za .- a) .+ mrl .^ 2)
end

function loop!(out, za, a, mrl)
    Nf, Nt = size(za)
    Na = length(a)
    Nm = length(mrl)
    @avx for m = 1:Nm
        for t = 1:Nt
            for j = 1:Na
                for i = 1:Nf
                    out[i, j, t, m] =
                        (1.0 - mrl[m]^2) /
                        (1.0 - 2.0 * mrl[m] * cos(za[i, t] - a[j]) + mrl[m]^2)
                end
            end
        end
    end
    return out
end

function compactloop!(out, za, a, mrl)
    Nf, Nt = size(za)
    Na = length(a)
    Nm = length(mrl)
    @avx for t = 1:Nt, j = 1:Na, i = 1:Nf, m = 1:Nm
        out[i, j, t, m] =
            (1.0 - mrl[m]^2) /
            (1.0 - 2.0 * mrl[m] * cos(za[i, t] - a[j]) + mrl[m]^2)
    end
    return out
end

function unroll!(out, za, a, mrl)
    Nf, Nt = size(za)
    Na = length(a)
    Nm = length(mrl)
    @avx for t = 1:Nt, j = 1:Na, i = 1:Nf
        c = cos(za[i, t] - a[j])
        out[i, j, t, 1] =
            (1.0 - mrl[1]^2) /
            (1.0 - 2.0 * mrl[1] * c + mrl[1]^2)
        out[i, j, t, 2] =
            (1.0 - mrl[2]^2) /
            (1.0 - 2.0 * mrl[2] * c + mrl[2]^2)
    end
    return out
end

function helper!(out, za, a)
    Nf, Nt = size(za)
    Na = length(a)
    @avx for t = 1:Nt, j = 1:Na, i = 1:Nf, m = 1:Nm
        out[i, j, t, m] = cos(za[i, t] - a[j])
    end
    return out
end

function base!(out)
    @avx for i in eachindex(out)
        out[i] = cos(out[i])
    end
end

function base2!(out)
    out .= cos.(out)
end


Nf = 266
Nt = 100
Na = 361
Nm = 2

za = 2π .* rand(Nf, Nt)
a = collect(range(0, 2π, length = Na))
mrl = rand(Nm)
out = Array{Float64}(undef, Nf, Na, Nt, Nm)
out2 = Array{Float64}(undef, Nf, Na, Nt)

@assert wrapcauchy!(out, za, a, mrl) == wrapcauchy(za, a, mrl)
@assert all(loop!(out, za, a, mrl) .≈ wrapcauchy(za, a, mrl))
@assert all(compactloop!(out, za, a, mrl) .≈ wrapcauchy(za, a, mrl))

println()
println("reference")
@btime out = wrapcauchy(za, a, mrl)
println("allocation")
@btime out = Array{Float64}(undef, Nf, Na, Nt, Nm)
println("inplace")
@btime wrapcauchy!(out, za, a, mrl)
println("loop")
@btime loop!(out, za, a, mrl)
println("compactloop")
@btime compactloop!(out, za, a, mrl)
println("unroll")
@btime compactloop!(out, za, a, mrl)
println("helper")
@btime helper!(out, za, a)
println("base")
@btime base!(out)
println("base2")
@btime base2!(out)
