using StructArrays
using LoopVectorization

function reference!(particles, i, j)
    @views particles[:, i] .= particles[:, j]
end

function unroll!(particles, i, j)
    @views particles.m[:, i] .= particles.m[:, j]
    @views particles.f[:, i] .= particles.f[:, j]
    @views particles.x[:, i] .= particles.x[:, j]
    @views particles.y[:, i] .= particles.y[:, j]
    @views particles.vx[:, i] .= particles.vx[:, j]
    @views particles.vy[:, i] .= particles.vy[:, j]
end

function replpacefieldscolumn!(particles, i, j)
    replacecolumn!(particles.m, i, j)
    replacecolumn!(particles.f, i, j)
    replacecolumn!(particles.x, i, j)
    replacecolumn!(particles.y, i, j)
    replacecolumn!(particles.vx, i, j)
    replacecolumn!(particles.vy, i, j)

end

function replacecolumn!(A, i, j)
    @inbounds for k = 1:size(A, 1)
        A[k, i] = A[k, j]
    end
end

Nt = 10000
Np = 10000
particles = StructArray{State}(undef, Nt, Np)
A = Array{Float64}(undef, Nt, Np)


println()
println("reference")
@btime reference!(particles, 1, Np)
println("unroll")
@btime unroll!(particles, 1, Np)
println("unrollfill")
@btime unrollfill!(particles, 1, Np)
