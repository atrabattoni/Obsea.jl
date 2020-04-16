function argsample(w, N)
    out = Array{Int64,1}(undef, N)
    j = 1
    s = w[j]
    u = rand() / N
    for i=1:N
        while s<u
            j += 1
            s += w[j]
        end
        out[i] = j
        u += 1 / N
    end
    out
end
