function likelihood(zr, tdoalut, models, grid)
    # parameters
    Nτ, Nt = size(zr)
    @unpack Nmode, v, τ = tdoalut
    @unpack Nr, Nm = grid
    @unpack lam = models
    # allocation
    ℓr = ones(Nr, Nt, Nm)
    u = Array{Float64}(undef, Nτ, Nt)
    out = Array{Float64}(undef, Nτ)
    for m = 1:Nm
        # non centered chi square
        u .= exp.(-lam[m] ./ 2.0) .* besseli.(0, sqrt.(lam[m] .* zr))
        for t = 1:Nt
            for mode = 1:Nmode
                # propagation incertainty
                @views convsame!(out, u[:, t], v[mode])
                # interpolation on the range grid
                itp = interpolate!(out, BSpline(Linear()))
                itp = extrapolate(itp, 1.0)
                itp = scale(itp, grid2range(grid.τ))
                ℓr[:, t, m] .*= itp(τ[:, mode])
            end
        end
    end
    return ℓr
end
function likelihood(za, models, grid)
    Nf, Nt = size(za)
    @unpack Na, Nm = grid
    @unpack mrl, n = models
    # wrapped cauchy distribution
    za = reshape(za, (Nf, 1, Nt, 1))
    a = reshape(grid.a, (1, Na, 1, 1))
    mrl = reshape(mrl, (1, 1, 1, Nm))
    ℓa = (1.0 .- mrl .^ 2) ./ (1.0 .- 2.0 .* mrl .* cos.(za .- a) .+ mrl .^ 2)
    # frequency stability
    out = Array{Float64}(undef, Nf)
    for m = 1:Nm
        for t = 1:Nt
            for j = 1:Na
                @views rollprod!(out, ℓa[:, j, t, m], n[m])
                ℓa[:, j, t, m] .= out
            end
        end
    end
    return ℓa
end

function marginalize(ℓr, ℓa, models, grid)
    # parameters
    @unpack Nr, Nf, Na, Nm = grid
    @unpack pb = models
    # marginalize
    pb = reshape(pb, (1, Nm))
    ℓm =
        pb .* dropdims(sum(grid.r .* ℓr, dims = 1), dims = 1) ./ sum(grid.r) .*
        dropdims(sum(ℓa, dims = (1, 2)), dims = (1, 2)) ./ Na ./ Nf
    ℓ0 = 1.0 - sum(pb)
    ℓΣm = ℓ0 .+ dropdims(sum(ℓm, dims = 2), dims = 2)
    return ℓm, ℓΣm
end

struct Likelihood
    r::Array{Float64,2}
    a::Array{Float64,3}
    m::Vector{Float64}
    Σm::Float64
end
