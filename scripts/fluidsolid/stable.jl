using BenchmarkTools

function Rpp(θ, ρ₁, ρ₂, α₁, α₂, β₂)
    x = sin(θ)
    if x > α₁ / α₂
        return 1.0
    end
    a = α₂ / α₁
    b = β₂ / α₂
    r = ρ₂ / ρ₁
    f(j, x) = sqrt(1 - (j * x)^2)
    g(j, x) = 1 - 2 * (j * x)^2
    c = -a * f(1, x) * g(b, x)^2 * r -
        4 * b^3 * f(1, x) * f(a, x) * f(b, x) * r * x^2
    d = 2 * b^2 * f(a, x) * x^2 + f(a, x) * g(b, x)
    return (c + d) / (c - d)
end


function propagation(r, z, ρ₁, ρ₂, α₁, α₂, β₂, n)
    k = range(1, 2 * n - 1, step = 2)
    k = collect(reshape(k, (1, length(k))))
    l = sqrt.((k .* z) .^ 2 .+ r .^ 2)
    t = l ./ α₁
    τ = diff(t, dims = ndims(t))
    θ = asin.(r ./ l)
    rpp = Rpp.(θ, ρ₁, ρ₂, α₁, α₂, β₂)
    a = rpp .^ k ./ l .^ 2
    a² = a .^2
    g = 2 * a[2:end] ./ a[1:end-1]
    @show θ rpp
    return τ, g
end

τ, g = propagation.(10000, 5000, 1000, 2000, 1500, 1800, 500, 7)
