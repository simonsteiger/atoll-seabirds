module CustomUtilityFuns

export showall,
       extractparams,
       pct,
       safelogistic,
       idx,
       sdim,
       rdistearth

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

function rdistearth(x1; x2=nothing, miles=false, R=nothing)
    isnothing(R) && begin R = miles ? 3963.34 : 6378.388 end
    coslat1, sinlat1 = [@. f((x1[:, 2] * π) / 180) for f in [cos, sin]]
    coslon1, sinlon1 = [@. f((x1[:, 1] * π) / 180) for f in [cos, sin]]
    
    if isnothing(x2)
        pp = [coslat1 .* coslon1 coslat1 .* sinlon1 sinlat1] * [coslat1 .* coslon1 coslat1 .* sinlon1 sinlat1]'
        return @. R * acos([abs(x > 1 ? 1 * sign(x) : x for x in pp)])
    else
        coslat2, sinlat2 = [@. f((x2[:, 2] * π) / 180) for f in [cos, sin]]
        coslon2, sinlon2 = [@. f((x2[:, 1] * π) / 180) for f in [cos, sin]]
        pp = [coslat1 .* coslon1 coslat1 .* sinlon1 sinlat1] * [coslat2 .* coslon2 coslat2 .* sinlon2 sinlat2]'
        return @. R * acos([abs(x > 1 ? 1 * sign(x) : x for x in pp)])
    end
end

end
