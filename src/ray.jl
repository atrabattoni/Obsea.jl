function tlim(incidence, depth, celerity, nmode)
    k = collect(1:2:2*nmode+1)
    t = k .* depth ./ cos(incidence) ./ celerity
end


function tdoa(r, depth, celerity, nmode, tc, tb)
    k = collect(1:2:2*nmode+1)
    t = zeros(nmode + 1)
    for i = 1:nmode+1
        t[i] = sqrt((k[i] * depth)^2 + r^2) / celerity
        if !(tc[i] < t[i] < tb[i])
            t[i] = NaN
        end
    end
    tau = diff(t)
end


function tau2range(depth, celerity, nmode, fs, nfft)

    dtau = 1.0 / fs
    ntau = nfft // 2 + 1
    tau = dtau * collect(1:ntau)

    f(tau, k) =
        sqrt(
            (tau * celerity - 2 * depth) *
            (tau * celerity + 2 * depth) *
            (tau * celerity - 4 * depth * k - 4 * depth) *
            (tau * celerity + 4 * depth * k + 4 * depth),
        ) / (2 * tau * celerity)

    r = zeros(nmode, ntau)
    for k = 1:nmode
        r[k] = f(tau, k)
    end

    return tau, r
end
