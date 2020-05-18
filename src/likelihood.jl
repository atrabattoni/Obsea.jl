"""
    likelihood(z, [tdoalut,] models, grid)

Compute the likelihood ratio on each point of grid for each models.

If tdoalut is given, the measurement z is supposed to be a ceptsrogram, and an
azigram otherwise.
"""
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
    ℓa = Array{Float64}(undef, Nf, Na, Nt, Nm)
    wrapcauchy!(ℓa, za, grid.a, mrl)
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

"""
    wrapcauchy!(out, za, a, mrl)

Compute the wrapcauchy probability divided by the uniform one.

Computation are done for each combination of measured angles za, hypothetic
angles a, and give mean running length mrl.
"""
function wrapcauchy!(out, za, a, mrl)
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

"""
    marginalize(ℓr, ℓa, models, grid)

Marginalize the likelihood to get the probability of each model.
"""
function marginalize(ℓr, ℓa, models, grid)
    # parameters
    @unpack Nr, Nf, Na, Nm = grid
    @unpack pb = models
    # marginalize
    ℓm =
        reshape(pb, (1, Nm)) .*
        dropdims(sum(grid.r .* ℓr, dims = 1), dims = 1) ./ sum(grid.r) .*
        dropdims(sum(ℓa, dims = (1, 2)), dims = (1, 2)) ./ Na ./ Nf
    ℓ0 = 1.0 - sum(pb)
    ℓΣm = ℓ0 .+ dropdims(sum(ℓm, dims = 2), dims = 2)
    return ℓm, ℓΣm
end

"""
    Likelihood(ℓr, ℓa, ℓm, ℓΣm)

Store all likelihood computations.
"""
struct Likelihood
    r::Array{Float64,3}
    a::Array{Float64,4}
    m::Array{Float64,2}
    Σm::Vector{Float64}
end
function Likelihood(zr, za, tdoalut, models, grid)
    ℓr = likelihood(zr, tdoalut, models, grid)
    ℓa = likelihood(za, models, grid)
    ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)
    return Likelihood(ℓr, ℓa, ℓm, ℓΣm)
end

"""
    LikelihoodSlice(ℓtr, ℓta, ℓtm, ℓtΣm)

Store all likelihood computations at a given time.
"""
struct LikelihoodSlice
    r::Array{Float64,2}
    a::Array{Float64,3}
    m::Vector{Float64}
    Σm::Float64
end
