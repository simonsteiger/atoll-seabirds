module Postprocess

export extractparams, prediction, targets

using Turing, StatsFuns

function pred_presence(p::NamedTuple, s, s_in_n, n, r, a, pc)
    y = fill(false, length(a))

    for i in eachindex(a)
        v = logistic(
            p.λ[r[i], s] +
            p.θ1[n][s_in_n] * pc[i, 1] +
            p.θ2[n][s_in_n] * pc[i, 2] +
            p.θ3[n][s_in_n] * pc[i, 3] +
            p.θ4[n][s_in_n] * pc[i, 4] +
            p.θ5[n][s_in_n] * pc[i, 5] +
            p.θ6[n][s_in_n] * pc[i, 6]
        )
        y[i] = rand(Bernoulli(v))
    end
    
    return y
end

function prediction(target, p::NamedTuple, s, s_in_n, n, r, a, pc)
    target == "presence" ? pred_presence(p::NamedTuple, s, s_in_n, n, r, a, pc) :
    target == "count" ? nothing :
    throw("Unknown target: $target")
end

end
