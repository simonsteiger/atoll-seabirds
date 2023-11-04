module CustomUtilityFuns

export showall,
       extractparams,
       pct,
       safelogistic,
       idx,
       sdim

using Turing, DataFrames

# Default method for printing all entries of an object
function showall(x)
    show(stdout, "text/plain", x)
end

# Method for printing all entries of the parameter panel of chains
function showall(x::Chains; table="parameters")
    idx = table == "parameters" ? 1 : 2
    x = describe(x)[idx]
    show(DataFrame(x), allrows=true)
end

# `extractparams()` is a wrapper around `get_params()` to allow easier iteration
extractparams(chain, param) = getproperty(get_params(chain), Symbol(param))

# `extractparams()` automatically loops over all entries of params if params is a Vector of Symbols or Strings
function extractparams(chain, params::AbstractArray{Symbol})
    values = [getproperty(get_params(chain), p) for p in params]
    namedvalues = (; zip(params, values)...)
    return namedvalues
end

function extractparams(chain, params::AbstractArray{String})
    sym_params = Symbol.(params)
    values = [getproperty(get_params(chain), p) for p in sym_params]
    namedvalues = (; zip(sym_params, values)...)
    return namedvalues
end

# Helper to print info message about divergent transitions
function pct(vec, val)
    round(sum(vec .== val) / length(vec) * 100, digits=2)
end

# Convert matrix indexing to vector indexing
idx(i, j) = i .+ (j .- 1) * maximum(i)

# This version of logistic should not underflow
safelogistic(x::T) where {T} = logistic(x) * (1 - 2 * eps(T)) + eps(T)

# Slice a dimension?
sdim(n) = (a) -> map(d -> d[n], a)

end
