module CustomUtilityFuns

using Turing, DataFrames

export getsamples,
       peaceful_generated_quantities,
       idx,
       lu,
       between,
       standardise, unstandardise,
       ModelSummary

# Shorthand length unique
lu(x) = length(unique(x))

# Convert matrix indexing to vector indexing
idx(i, j) = i .+ (j .- 1) * maximum(i)

# Z standardise
standardise(v) = (v .- mean(v)) ./ std(v)
# Revert Z standardisation
unstandardise(z, v) = z .* std(v) .+ mean(v)

# Get samples from a generated quantities object
getsamples(θ, fields...) = map(x -> NamedTuple{fields}(getfield.(Ref(x), fields)), θ)
# TODO make a generated quantities struct?

# Get generated quantities from model / chain and avoid warning for model internals
peaceful_generated_quantities(m, c) = generated_quantities(m, Turing.MCMCChains.get_sections(c, :parameters))

# Summarise a DynamicPPL Model
struct ModelSummary{T}
    model::T
    chains::Chains
    θ
    samples
    function ModelSummary(model, chains)
        quantities = peaceful_generated_quantities(model, chains)
        θ = keys(quantities[1])
        samples = vec(getsamples(quantities, θ...))
        return new{typeof(model)}(model, chains, θ, samples)
    end
end

between(x, lower, upper) = lower ≤ x && x ≤ upper

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
