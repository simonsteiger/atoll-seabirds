# This script is part of the project associated with the article
# ...
# Authors:
# Last edited:

"This module contains helper functions for modeling and postprocessing."
module CustomUtilityFuns

# Modeling and wrangling
using Turing, DataFrames, Chain
# Plotting
using Colors, ColorSchemes, StatsPlots

export getsamples,
       peaceful_generated_quantities,
       idx, midx,
       lu,
       between,
       standardise, unstandardise,
       ModelSummary,
       diagnose,
       popsum,
       calcratio,
       blue_in_blue

"Shorthand length unique"
lu(x) = length(unique(x))

"Convert matrix indexing to vector indexing"
idx(i, j) = i .+ (j .- 1) * maximum(i)

"Convert vector indexing to matrix indexing"
function midx(x, r)
    i = x % r == 0 ? r : x % r
    j = i == r ? x ÷ r : x ÷ r + 1
    return [j, i]
end

"Z standardise"
standardise(v) = (v .- mean(v)) ./ std(v)
"Revert Z standardisation"
unstandardise(z, v) = z .* std(v) .+ mean(v)

"Get samples from a generated quantities object"
getsamples(θ, fields...) = map(x -> NamedTuple{fields}(getfield.(Ref(x), fields)), θ)

"Get generated quantities from model / chain and avoid warning for model internals"
peaceful_generated_quantities(m, c) = generated_quantities(m, Turing.MCMCChains.get_sections(c, :parameters))

"Summarise a DynamicPPL Model"
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

"Calculate ratio from an array of binary values"
ratio(x::AbstractArray) = sum(x) / length(x)

"Print summary of basic MCMC diagnostics"
function diagnose(chains::Chains)
    numerical_error_ratio = ratio(DataFrame(chains).numerical_error)
    df_describe = DataFrame(describe(chains)[1])
    max_rhat = maximum(df_describe.rhat)
    min_ess_bulk, min_ess_tail = [minimum(col) for col in eachcol(df_describe[:, [:ess_bulk, :ess_tail]])]
    
    msg = """Minimal MCMC Diagnostics
    Divergent transitions    => $(round(numerical_error_ratio*100, digits=2))%
    Maximum Rhat             => $(round(max_rhat, digits=4))
    Minimum ESS (Bulk, Tail) => $(Int64(round(min_ess_bulk, digits=0))), $(Int64(round(min_ess_tail, digits=0)))
    """

    if numerical_error_ratio > 0.03 || max_rhat > 1.01 || any(.>(1000, [min_ess_bulk, min_ess_tail]))
        @warn msg
    else
        @info msg
    end
        
    return nothing
end

"Wrapper around simple summary for easier broadcasting in count summary"
function popsum(df)
    out = @chain df begin
        groupby(_, :species)
        combine(_, [:nbirds, :lower, :upper] .=> sum .=> identity)
    end
    return out
end

"Compare populations on atolls to global population of each species"
function calcratio(species, n, blmin, blmax, hbw, otero)
    species ∈ ["Gygis_alba", "Nesofregetta_fuliginosa"] ? n / hbw :
    ismissing(blmin) && ismissing(hbw) ? n / otero :
    ismissing(blmin) ? n / hbw :
    ismissing(blmax) ? n / blmin :
    n / mean([blmin, blmax])
end

"Return pseudo-categorical color scheme for two cutoffs"
function in_out_colorscheme(out, in; eps=eps(Float64))
    return cgrad([out, in, in, out], [0.0, eps, 1.0 - eps, 1.0])
end

# Create color scheme for validation plots in count pipeline
out_blue = HSLA(200, 0.9, 0.8, 0.8)
in_blue = HSLA(200, 0.9, 0.4, 0.8)
blue_in_blue = in_out_colorscheme(out_blue, in_blue)

end
