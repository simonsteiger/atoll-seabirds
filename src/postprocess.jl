module Postprocess

export extractparams, prediction, targets

using Turing

# `extractparams()` is a wrapper around `get_params()` to allow easier iteration
extractparams(chain, param) = getproperty(get_params(chain), Symbol(param))

# `extractparams()` automatically loops over all entries of params if params is a Vector of Symbols or Strings
function extractparams(chain, params::AbstractArray{T}) where {T<:Union{String,Symbol}}
    values = [getproperty(get_params(chain), p) for p in Symbol.(params)]
    namedvalues = (; zip(params, values)...)
    return namedvalues
end

omegas = [Symbol("ω$i$j") for i in 1:6, j in 1:3]
targets = [:α, :β, :γ, omegas...]

d = extractparams(chain, targets)
d.ω11

keys(d)

function pred_presence(p::NamedTuple, o::NamedTuple)
    presence = fill(zeros(length(o.atoll)))

    for i in eachindex(o.atoll)
        v = logistic(
            p.α[o.atoll[i]] +
            p.β[o.species[i]] +
            p.γ[o.atoll[i]] + # not sure if this enters the model like that
            p.ω1[o.nesting[i]][o.species[i]] * o.pc[i, 1] +
            p.ω2[o.nesting[i]][o.species[i]] * o.pc[i, 2] +
            p.ω3[o.nesting[i]][o.species[i]] * o.pc[i, 3] +
            p.ω4[o.nesting[i]][o.species[i]] * o.pc[i, 4] +
            p.ω5[o.nesting[i]][o.species[i]] * o.pc[i, 5] +
            p.ω6[o.nesting[i]][o.species[i]] * o.pc[i, 6]
        )
        presence[i] = Bernoulli(v)
    end

    return presence
end

function prediction(target, p::NamedTuple, o::NamedTuple)
    target == "presence" ? pred_presence(p, o) :
    target == "count"    ? nothing             :
    throw("Unknown target: $target")
end
    
end
