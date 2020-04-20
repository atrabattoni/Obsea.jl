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
