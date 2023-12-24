module PresenceModel

export nothing

# --- WORKSPACE SETUP --- #

# Pathing
const ROOT = dirname(Base.active_project())

# Probabilistic programming
using Turing, ReverseDiff, ParetoSmooth
# Benchmarking
using TuringBenchmarking, BenchmarkTools
# Model speed optimization
using LazyArrays
# Statistics
using StatsFuns, LinearAlgebra
# Working with tabular data
using Chain, DataFrames
# Plotting
using StatsPlots
# Saving results and logging
using Serialization, CSV, Dates, Markdown, OrderedCollections
# Random seeds
using Random

# Load custom modules
include("$ROOT/src/presence.jl")
using .PresenceVariables
include("$ROOT/src/utilities.jl")
using .CustomUtilityFuns
include("$ROOT/scripts/modelzoo.jl")
using .PresenceModels

# Set seed
Random.seed!(42)

# --- MODEL SETUP --- #

model = mMp
simulate = simulate_Mp
odict_inputs = OrderedDict(
    "region" => region.known.num,
    "species" => species.known.num,
    "nesting" => nesting.known.num,
    "PC" => PC.known,
    "species_within_nesting" => species_in_nesting.known,
    "unique_nesting" => nesting.levels,
    "unique_species_within_nesting" => species_in_nesting.levels,
    "n_burrow" => nspecies.burrow,
    "n_ground" => nspecies.ground,
    "n_vegetation" => nspecies.vegetation,
)

# Path to read or write chain from / to
chainpath = "presence.jls"

# --- PRIOR PREDICTIVE CHECK --- #

m = model(values(odict_inputs)..., presence)

prior = @chain begin
    sample(m, Prior(), 1000)
    ModelSummary(m, _)
end

let
    priorpreds = reduce(vcat, simulate(prior.samples, values(odict_inputs)...))
    density(priorpreds, normalize=true, c=2, fillrange=0, fillalpha=0.2, legend=false)
end

# --- MODEL CONFIG --- #

backends = [
    Turing.Essential.ForwardDiffAD{0}(),
    Turing.Essential.ReverseDiffAD{false}(),
    Turing.Essential.ReverseDiffAD{true}()
];

TuringBenchmarking.run(TuringBenchmarking.make_turing_suite(m, adbackends=backends);)
@info "TuringBenchmark shows that ReverseDiff{true} is the fastest AD backend."

# Use ReverseDiffAD and enable caching
Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

# Configure sampling
sampler = NUTS(1000, 0.95; max_depth=10)
nsamples = 2000
nchains = 4
config = (sampler, MCMCThreads(), nsamples, nchains)

@info """Sampler: $(string(sampler))
Samples: $(nsamples)
Threads: $(nchains)
"""

@info "🚀 Starting sampling: $(Dates.now())"
posterior = @chain begin
    sample(m, config...)
    ModelSummary(m, _)
end

try
    path = "$ROOT/results/chains/$chainpath"
    serialize(path, posterior.chains)
    @info "Saved chains to `$path`."
catch error
    @warn "Writing chains failed with an $error."
end

# --- POSTERIOR PREDICTIVE CHECK --- #

postpcplots = let preds = reduce(hcat, simulate(posterior.samples, values(odict_inputs)...))
    # Iterate over species and create plots with posterior predictions
    map(enumerate(unique(species.known.num))) do (idx, sp)
        # Predicted values for species sp
        pred_x = vec(preds[species.known.num.==sp, :])
        # Observed values for species sp
        obs_x = presence[species.known.num.==sp]
        # Assemble plot for species sp
        histogram(obs_x, normalize=true, alpha=0.5, lc=:transparent, bins=10, label="O")
        density!(pred_x, fillrange=0, fillalpha=0.2, normalize=true, alpha=0.5, bins=10, lc=:transparent, yticks=:none, label="P")
        xticks!(0:0.5:1, string.(0:0.5:1))
        title!(unique(species.known.str)[idx], titlefontsize=8)
        vline!([0.8], c=:black, ls=:dash, label=:none)
    end
end

plot(postpcplots..., titlefontsize=9, size=(800, 1200), layout=(8, 5))
savefig("results/svg/postpreds_Mp.svg")

# Create predictions for unknown atolls
predictions = let odict_unknown_inputs = copy(odict_inputs)
    names = ["region", "species", "PC", "species_within_nesting"]
    vals = [region.unknown.num, species.unknown.num, PC.unknown, species_in_nesting.unknown]
    setindex!.(Ref(odict_unknown_inputs), vals, names)
    preds = simulate(posterior.samples, values(odict_unknown_inputs)...)
    reduce(hcat, preds)
end

# Wrap predictions and quantile in DataFrame with inputs
df_preds = let
    μs = vec(mean(predictions, dims=2))
    qs = reduce(hcat, [quantile(slice, [0.05, 0.95]) for slice in eachslice(predictions, dims=1)])
    DataFrame(
        [atoll.unknown.str, region.unknown.str, species.unknown.str, μs, qs[1, :], qs[2, :]],
        [:atoll, :region, :species, :mean, :lower005, :upper095]
    )
end


# Plot predictions for unknown atolls and store results in Dict
dict_unknown_pred_plots = let df = df_preds
    plots = map(enumerate(unique(df.species))) do (index, species)
        I = df.species .== species
        scatter(df.mean[I], df.atoll[I], ms=2, msc=1, xerror=(abs.(df.lower005[I] .- df.mean[I]), abs.(df.upper095[I] .- df.mean[I])))
        title!(unique(species.unknown.str)[index], titlefontsize=8)
        xlims!(0, 1)
        vline!([0.8], c=:black, ls=:dash, legend=:none)
    end
    Dict(Pair.(unique(species.unknown.str), plots))
end

plot(dict_unknown_pred_plots..., titlefontsize=9, size=(800, 1200), layout=(8, 5))
# TODO plot by nesting type or genus

# Save predictions for unknown atolls to CSV
try
    path = "$ROOT/results/data/presencepreds.csv"
    CSV.write(path, df_preds); @info "Saved predictions to `$path`."
catch error
    @warn "Writing predictions failed with an $error."
end

# LOO
cv_res = psis_loo(m, posterior.chains);
cv_res
# - no overfit
# - out of sample performance near in-sample performance (gmpd 0.76)
# - not many outliers in the p_eff plot (the outliers are logical => Clipperton, Ant (PCs?))
# - in line with posterior predictive check
# - not sure how to interpret differences in naive_lpd and cv_elpd ... but they seem very low p_avg 0.02

end
