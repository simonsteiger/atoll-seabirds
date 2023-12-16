module CustomUtilityFuns

export showall,
       getsamples,
       peaceful_generated_quantities,
       idx,
       lu,
       rdistearth

using Turing, DataFrames

# TODO split this into modeling related utilities like idx, lu (not included yet), getsamples...
# and others like ... ? rdistearth, showall? But showall isn't really very important in the final version :thinking:

# Default method for printing all entries of an object
function showall(x)
    show(stdout, "text/plain", x)
end

lu(x) = length(unique(x))

# Method for printing all entries of the parameter panel of chains
function showall(x::Chains; table="parameters")
    idx = table == "parameters" ? 1 : 2
    x = describe(x)[idx]
    show(DataFrame(x), allrows=true)
end

function peaceful_generated_quantities(m, c)
    chains_params = Turing.MCMCChains.get_sections(c, :parameters)
    return generated_quantities(m, chains_params)
end

getsamples(θ, fields...) = map(x -> NamedTuple{fields}(getfield.(Ref(x), fields)), θ)

# Convert matrix indexing to vector indexing
idx(i, j) = i .+ (j .- 1) * maximum(i)

# This version of logistic should not underflow
safelogistic(x::T) where {T} = logistic(x) * (1 - 2 * eps(T)) + eps(T)

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
