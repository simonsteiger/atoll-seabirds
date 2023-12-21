# --- WORKSPACE SETUP --- #

# Probabilistic programming
using Turing, TuringBenchmarking, ReverseDiff, ParetoSmooth
# Model speed optimization
using LazyArrays
# Statistics
using StatsFuns
# Covariance matrices
import LinearAlgebra: I
# Storing and working with data
using Chain, DataFrames, OrderedCollections
# Plotting
using StatsPlots
# Saving results and logging
using Serialization, CSV, Dates, Markdown
# Random seeds
using Random

# Path
const ROOT = dirname(Base.active_project())

# Load custom modules
include("$ROOT/scripts/modelzoo.jl")
include("$ROOT/src/count.jl")
include("$ROOT/src/utilities.jl")

# Make custom modules available
using .CountVariables
using .CountModels
using .CustomUtilityFuns

# Benchmark model?
benchmark = false

# Run the sampler?
run = isempty(ARGS) ? true : ARGS[1] == "true"

# Pick model from ModelZoo 
# We might as well paste the model here at the end and just reference the model zoo for what we tried?
model = mMc
simulate = simulate_Mc
# Set inputs
odict_inputs = OrderedDict(
    "region" => num_region_known,
    "species" => num_species_known,
    "nesting" => num_nesting_known,
    "PC" => PC_known,
    "species_within_nesting" => num_species_within_nesting_known,
    "unique_nesting" => unique_nesting_known,
    "unique_species_within_nesting" => unique_species_within_nesting_known,
    "n_burow" => count_species_by_nesting[1],
    "n_ground" => count_species_by_nesting[2],
    "n_vegetation" => count_species_by_nesting[3],
)

# If not loading a chain, save results to path below
chainpath = "count.jls"

zlogn = standardise(log.(nbirds))

# --- PRIOR PREDICTIVE CHECKS --- #

m = model(values(odict_inputs)..., zlogn)

prior = @chain begin
    sample(m, Prior(), 1000)
    ModelSummary(m, _)
end

let
    priorpreds = reduce(vcat, simulate(prior.samples, values(odict_inputs)...))
    histogram(zlogn, normalize=true, c=1)
    density!(priorpreds, normalize=true, c=:white, lw=3)
    density!(priorpreds, normalize=true, c=2, fillrange=0, fillalpha=0.2, legend=false)
    xlims!(-5, 5)
end
# not sure if necessary to show both observed and prior for prior predictive check
# I guess it's not bad on a second thought though – we don't want probability mass in useless places

# Benchmark different backends to find out which is fastest

backends = [
    Turing.Essential.ForwardDiffAD{0}(),
    Turing.Essential.ReverseDiffAD{false}(),
    Turing.Essential.ReverseDiffAD{true}()
]

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

# Notify user about sampling configuration
@info """Sampler: $(string(sampler))
Samples: $(nsamples)
Chains: $(nchains)
"""

# Set seed
Random.seed!(42)

# Discarding samples is unnecessary after NUTS tuning
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

posteriorpreds = simulate(posterior.samples, values(odict_inputs)...)

let preds = reduce(vcat, posteriorpreds)
    histogram(zlogn, normalize=true, c=1)
    density!(preds, normalize=true, c=:white, lw=3)
    density!(preds, normalize=true, c=2, fillrange=0, fillalpha=0.2, legend=false)
    xlims!(-5, 5)
end

# Create dictionary of species-wise posterior prediction plot
dict_postpred_plot = let
    # Index vector I sorts other vectors to match region order
    I = sortperm(num_region_known)
    # Create matrix of posterior predictions, rows X is species X
    # Unstandardised and on count scale
    preds = exp.(unstandardise(reduce(hcat, posteriorpreds), log.(nbirds)))
    # Enumerate species keys to create one plot per species
    plots = map(enumerate(keys(odict_species))) do (index, species)
        # Index for subsetting to species S in current iteration
        S = num_species_known[I] .== index
        # Get observed values for species X
        obs = nbirds[I][S]
        # Calculate mean of predicted values for species X
        pred = vec(mean(preds[I, :][S, :], dims=2))
        # Calculate quantile of predicted values for species X
        qlims = @chain preds[I, :][S, :] begin
            map(x -> quantile(x, [0.025, 0.975]), eachslice(_, dims=1))
            map(i -> abs.(_[i] .- pred[i]), eachindex(_))
            reduce(hcat, _)
        end
        # Assemble plot for species X
        scatter(str_atoll_known[I][S], pred, yerror=(qlims[1, :], qlims[2, :]), markersize=2.5, msc=1, label="P", permute=(:y, :x))
        scatter!(str_atoll_known[I][S], obs, markersize=2.5, msc=2, label="O", permute=(:y, :x))
        title!(replace(species, "_" => " "), titlefontsize=12)
    end
    # Collect all plots in a dictionary (plotting all at once is a bit busy, maybe by nesting type or genus?)
    Dict(Pair.(keys(odict_species), plots))
end

# TODO what's up with A minutus? Can we tweak the model to make it work better?

predictions =
    let thresholds = 0.75:0.05:0.85
        odict_unknown_inputs = copy(odict_inputs)
        names = ["region", "species", "PC", "species_within_nesting"]
        preds = map(thresholds) do threshold
            pass = ppres .> threshold
            vals = [num_region_unknown[pass], num_species_unknown[pass], PC_unknown[pass, :], num_species_within_nesting_unknown[pass]]
            setindex!.(Ref(odict_unknown_inputs), vals, names)
            samples = reduce(hcat, simulate(posterior.samples, values(odict_unknown_inputs)...))
        end
        Dict(Pair.(thresholds, preds))
    end

# TODO create input dict, change input dict here with oos vals, go
countpreds_oos = @chain begin
    predictcount(
        α,
        β,
        σ2,
        num_species_within_nesting_oos,
        num_species_oos,
        num_region_oos,
        PC_oos,
    )
    reduce(hcat, _)
    Matrix{Float64}(_)
    mean(_, dims=2)
    exp.(_)
    #[between(_[i], oos_lims[i][1], oos_lims[i][2]) for i in eachindex(_)]
    [_ oos_lims [between(_[i], oos_lims[i][1], oos_lims[i][2]) for i in eachindex(_)] sort(unique(str_species_known))[num_species_oos]]
end

avg_preds_unknown = vec(mean(countpreds_unknown, dims=2))
avg_preds_known = vec(mean(countpreds_known, dims=2))

df_countpreds = DataFrame(
    [
        num_atoll_unknown[threshold],
        num_region_unknown[threshold],
        num_species_unknown[threshold],
        avg_preds_unknown
    ],
    [:atoll, :region, :species, :nbirds]
)

CSV.write("$ROOT/results/data/countpreds_$PRIORSUFFIX.csv", df_countpreds)
@info "Successfully saved predictions to `$ROOT/results/data/countpreds_$PRIORSUFFIX.csv`."

# --- LOO CV --- #

# It seems that specifying the model as MvNormal confuses psis_loo ("1 data point")
# Respecify the likelihood as vectorized Normals
# y .~ Normal.(μ, σ)

# Regarding the "let's just put the model back into this script" ... 
# Not a big issue at all, but spelling out loospecial again seems nasty

loomodel = loospecial(inputs..., zlogn);

cv_res = psis_loo(loomodel, chain)
