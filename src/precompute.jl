struct Likelihood
    r::Array{Float64,2}
    a::Array{Float64,3}
    m::Vector{Float64}
    Σm::Float64
end

function likelihood(zr, za, tdoalut, models, grid)
    ℓr = likelihood(zr, tdoalut, models, grid)
    ℓa = likelihood(za, models, grid)
    ℓm, ℓΣm = marginalize(ℓr, ℓa, models, grid)
    return Likelihood(ℓr, ℓa, ℓm, ℓΣm)
end
function likelihood(z, tdoalut, models, grid)
    @unpack Nmode, v, τ = tdoalut
    @unpack Nr, Nm = grid
    ℓ = ones(Nr, Nm)
    @inbounds for m = 1:Nm
        @unpack lam = models[m]
        u = exp.(-lam ./ 2.0) .* besseli.(0, sqrt.(lam .* z))
        for mode = 1:Nmode
            A = convsame(u, v[mode])
            itp = interpolate(A, BSpline(Linear())) #
            itp = extrapolate(itp, 1.0)
            itp = scale(itp, grid2range(grid.τ))
            ℓ[:, m] .*= itp.(τ[:, mode])
        end
    end
    return ℓ
end
function likelihood(z, models, grid)
    @unpack Nf, Na, Nm = grid
    ℓ = zeros(Nf, Na, Nm)
    @inbounds for m = 1:Nm
        @unpack mrl, n = models[m]
        for (j, a) in enumerate(grid.a)
            for i = 1:Nf
                ℓ[i, j, m] = wrapcauchy(z[i], a, mrl)
            end
            ℓ[:, j, m] = rollprod(ℓ[:, j, m], n)
        end
    end
    return ℓ
end

function marginalize(ℓr, ℓa, models, grid)
    @unpack Nr, Nf, Na, Nm = grid
    pb = [model.pb for model in models]
    ℓm = similar(pb)
    @inbounds for m = 1:Nm
        @views ℓm[m] =
            pb[m] * sum(grid.r .* ℓr[:, m]) / sum(grid.r) * sum(ℓa[:, :, m]) /
            Na / Nf
    end
    ℓ0 = 1.0 - sum(pb)
    ℓΣm = ℓ0 + sum(ℓm)
    return ℓm, ℓΣm
end
