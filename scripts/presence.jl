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
include("$ROOT/src/utilities.jl")
include("$ROOT/scripts/modelzoo.jl")

# Make custom modules available
using .PresenceVariables
using .CustomUtilityFuns
using .PresenceModels

# Set seed
Random.seed!(42)

# --- SETTINGS --- #

model = mMp
simulate = simulate_Mp
odict_inputs = OrderedDict(
    "region" => num_region,
    "species" => num_species,
    "nesting" => num_nesting,
    "PC" => PC,
    "species_within_nesting" => num_species_within_nesting,
    "unique_nesting" => unique_nesting,
    "unique_species_within_nesting" => unique_species_within_nesting,
    "n_burow" => count_species_by_nesting[1],
    "n_ground" => count_species_by_nesting[2],
    "n_vegetation" => count_species_by_nesting[3],
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

backends = [
    Turing.Essential.ForwardDiffAD{0}(),
    Turing.Essential.ReverseDiffAD{false}(),
    Turing.Essential.ReverseDiffAD{true}()
];

TuringBenchmarking.run(TuringBenchmarking.make_turing_suite(m, adbackends=backends);)
@info "ReverseDiff{true} is the fastest AD backend."

# Set AD backend to :reversediff and compile with setrdcache(true)
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
    serialize("$ROOT/results/chains/$chainpath", posterior.chains); @info "💾 Chain saved to `$ROOT/results/chains/$chainpath`."
catch error
    @warn "Writing failed with an $error."
end

# --- POSTERIOR PREDICTIVE CHECK --- #

postpcplots = let
    postsamples = reduce(hcat, simulate(posterior.samples, values(odict_inputs)...))

    map(enumerate(unique(num_species))) do (index, species)
        pred_x = vec(postsamples[num_species.==species, :])
        obs_x = presence[num_species.==species]
        histogram(obs_x, normalize=true, alpha=0.5, lc=:transparent, bins=10, label="O")
        density!(pred_x, fillrange=0, fillalpha=0.2, normalize=true, alpha=0.5, bins=10, lc=:transparent, yticks=:none, label="P")
        xticks!(0:0.5:1, string.(0:0.5:1))
        title!(unique(str_species)[index], titlefontsize=8)
        vline!([0.8], c=:black, ls=:dash, label=:none)
    end
end

plot(postpcplots..., titlefontsize=9, size=(800, 1200), layout=(8, 5))
savefig("results/svg/postpreds_Mp.svg")

# Expand vectors to correct lenght
nau_long = [fill.(num_atoll_unknown, length(num_nesting_unknown))...;]
nru_long = [fill.(num_region_unknown, length(num_nesting_unknown))...;]
nsu_long = [fill(num_species_unknown, length(num_region_unknown))...;]
nnu_long = [fill(num_species_within_nesting_unknown, length(num_region_unknown))...;]
PCu_long = reduce(vcat, [permutedims(hcat(fill(s, length(num_nesting_unknown))...)) for s in eachslice(PC_unknown, dims=1)])

predictions = let odict_unknown_inputs = copy(odict_inputs)
    names = ["region", "species", "PC", "species_within_nesting"]
    vals = [nru_long, nsu_long, PCu_long, nnu_long]
    setindex!.(Ref(odict_unknown_inputs), vals, names)
    preds = simulate(posterior.samples, values(odict_unknown_inputs)...)
    reduce(hcat, preds)
end

df_preds = let 
    μs = vec(mean(predictions, dims=2))
    qs = reduce(hcat, [quantile(slice, [0.05, 0.95]) for slice in eachslice(predictions, dims=1)])
    DataFrame([nau_long, nru_long, nsu_long, μs, qs[1, :], qs[2, :]], [:atoll, :region, :species, :mean, :lower005, :upper095])
end


unknown_pred_plots = let df = df_preds
    plots = map(enumerate(unique(df.species))) do (index, species)
        I = df.species.==species
        scatter(df.mean[I], df.atoll[I], ms=2, msc=1, xerror=(abs.(df.lower005[I].-df.mean[I]), abs.(df.upper095[I].-df.mean[I])))
        title!(unique(str_species)[index], titlefontsize=8)
        xlims!(0, 1); vline!([0.8], c=:black, ls=:dash, legend=:none)
    end
    Dict(Pair.(unique(str_species), plots))
end

plot(unknown_pred_plots..., titlefontsize=9, size=(800, 1200), layout=(8, 5))

CSV.write("$ROOT/results/data/presencepreds_$PRIORSUFFIX.csv", df_preds)
@info "Successfully saved predictions to `$ROOT/data/presencepreds_$PRIORSUFFIX.csv`."

## Posterior predictive checks

posteriorpreds = Matrix{Float64}(reduce(hcat, vec(predictpresence(α, β, num_region, num_species, num_nesting, PC));))

# how do we best visualize binary posterior predictions?

scatter(mean.([presence, posteriorpreds]), yerror=std(presence))

# LOO
cv_res = psis_loo(model, chain);
println(cv_res)
# - no overfit
# - out of sample performance near in-sample performance (gmpd 0.76)
# - not many outliers in the p_eff plot (the outliers are logical => Clipperton, Ant (PCs?))
# - in line with posterior predictive check
# - not sure how to interpret differences in naive_lpd and cv_elpd ... but they seem very low p_avg 0.02

# TODO compare model L and model X - X looks better from the posterior predicitve plots, what does PSIS LOO say?