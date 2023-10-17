module Postprocess

export extractparams, prediction, targets

using Turing, StatsFuns

function pred_presence(p::NamedTuple, s, s_in_n, n, r, a, pc)
    y = fill([false, true], length(unique(a)))
    k = maximum(s)

    θ1n = getindex(p, Symbol("θ1$n"))[s_in_n]
    θ2n = getindex(p, Symbol("θ2$n"))[s_in_n]
    θ3n = getindex(p, Symbol("θ3$n"))[s_in_n]
    θ4n = getindex(p, Symbol("θ4$n"))[s_in_n]
    θ5n = getindex(p, Symbol("θ5$n"))[s_in_n]
    θ6n = getindex(p, Symbol("θ6$n"))[s_in_n]

    for i in eachindex(a)
        v = logistic.(
            p.λ[s+(r[i]-1)*k] + # Translate matrix indexing to vector indexing
            θ1n * pc[i, 1] +
            θ2n * pc[i, 2] +
            θ3n * pc[i, 3] +
            θ4n * pc[i, 4] +
            θ5n * pc[i, 5] +
            θ6n * pc[i, 6]
        )
        # v is a Matrix, but should be a vector
        y[i] = rand.(Bernoulli.(vec(v)))
    end

    return y
end

function prediction(target, p::NamedTuple, s, s_in_n, n, r, a, pc)
    target == "presence" ? pred_presence(p::NamedTuple, s, s_in_n, n, r, a, pc) :
    target == "count" ? nothing :
    throw("Unknown target: $target")
end

end
