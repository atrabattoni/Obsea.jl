function tlim(incidence, nmode, depth, celerity)
    k = collect(1:2:2*nmode+1)
    t = k .* depth ./ cos(incidence) ./ celerity
end


function tdoa(r, nmode, depth, celerity, tc, tb)
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


function convsame(u, v)
    @assert length(u) > length(v)
    @assert isodd(length(v))
    padding = length(v) ÷ 2
    conv(u, v)[1+padding:end-padding]
end

