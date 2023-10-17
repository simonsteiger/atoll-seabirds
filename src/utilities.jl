module CustomUtilityFuns

export showall,
       extractparams

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

lu(x) = length(unique(x))

end
