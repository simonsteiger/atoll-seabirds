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
include("$ROOT/src/count.jl")
using .CountVariables
include("$ROOT/scripts/modelzoo.jl")
using .CountModels
include("$ROOT/src/utilities.jl")
using .CustomUtilityFuns

# --- MODEL SETUP --- #

# Pick model from ModelZoo 
# We might as well paste the model here at the end and just reference the model zoo for what we tried?
model = mMc
simulate = simulate_Mc
# Set inputs
odict_inputs = OrderedDict(
    "region" => region.known.num,
    "species" => species.known.num,
    "nesting" => nesting.known.num,
    "PC" => PC.known,
    "species_within_nesting" => species_in_nesting.known.num,
    "unique_nesting" => nesting.levels,
    "unique_species_within_nesting" => species_in_nesting.levels,
    "n_burow" => nspecies.burrow,
    "n_ground" => nspecies.ground,
    "n_vegetation" => nspecies.vegetation,
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

# --- MODEL CONFIG --- #

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

# --- SAMPLING --- #

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

preds_train = simulate(posterior.samples, values(odict_inputs)...)

let preds = reduce(vcat, preds_train)
    histogram(zlogn, normalize=true, c=1)
    density!(preds, normalize=true, c=:white, lw=3)
    density!(preds, normalize=true, c=2, fillrange=0, fillalpha=0.2, legend=false)
    xlims!(-5, 5)
end

# Create dictionary of species-wise posterior prediction plot
dict_postpred_plot = 
    let preds = exp.(unstandardise(reduce(hcat, preds_train), log.(nbirds)))
        sorted_species = sort(unique(species.known.str))
        # I reorders other vectors by region
        I = sortperm(region.known.num)

        # Enumerate sorted species to create one plot per species
        plots = map(enumerate(sorted_species)) do (idx, sp)
            # Create BitVector S for species subsetting
            S = species.known.num[I] .== idx
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
            scatter(atoll.known.str[I][S], pred, yerror=(qlims[1, :], qlims[2, :]), markersize=2.5, msc=1, label="P", permute=(:y, :x))
            scatter!(atoll.known.str[I][S], obs, markersize=2.5, msc=2, label="O", permute=(:y, :x))
            title!(replace(sp, "_" => " "), titlefontsize=12)
        end
        # Collect all plots in a dictionary (plotting all at once is a bit busy, maybe by nesting type or genus?)
        Dict(Pair.(sorted_species, plots))
    end

# --- VALIDATION --- #

# The data below consists of islands where only population ranges are known in the literature.
# We make predictions for these data and check how much of the posterior lies within these population ranges.

preds_validation = let odict_oos_inputs = copy(odict_inputs)
    # Set names and values of inputs to be replaced in input dictionary
    names = ["region", "species", "PC", "species_within_nesting"]
    vals = [region.validation.num, species.validation.num, PC.validation, species_in_nesting.validation]
    setindex!.(Ref(odict_oos_inputs), vals, names)
    # Make predictions for validation data, transform to natural scale, and plot
    @chain posterior.samples begin
        simulate(_, values(odict_oos_inputs)...)
        reduce(hcat, _)
        exp.(unstandardise(_, log.(nbirds)))
        map(enumerate(eachslice(_, dims=1))) do (i, slice)
            histogram(slice, title=species.validation.str[i])
            xlims!(0, quantile(slice, 0.95))
            vline!(oos_lims[i], lw=2, c=:red)
        end
    end
end

# --- TARGET PREDICTIONS --- #

# The data below are those atolls where no information on count data are present.
# We predict data for these atolls, using three different thresholds, and summarise it for export.

preds_target =
    let thresholds = 0.75:0.05:0.85, odict_unknown_inputs = copy(odict_inputs)
        # Set names to be replaced in input dict, values will be subset in each iteration
        names = ["region", "species", "PC", "species_within_nesting"]
        preds = map(thresholds) do threshold
            # Create BitVector with presence predictions above threshold for indexing input vectors
            pass = ppres .> threshold
            vals = [region.unknown.num[pass], species.unknown.num[pass], PC.unknown[pass, :], species_in_nesting.unknown[pass]]
            setindex!.(Ref(odict_unknown_inputs), vals, names)
            # Predict values for those atolls above threshold, and transform to natural scale
            samples = @chain posterior.samples begin 
                simulate(_, values(odict_unknown_inputs)...)
                reduce(hcat, _)
            end
            # Calculate rowwise (atoll × species) mean and quantiles (I still haven't looked into why the correct dim here is 2)
            μs = exp.(unstandardise(mean(samples, dims=2), log.(nbirds)))
            qs = @chain samples begin
                [quantile(slice, [0.05, 0.95]) for slice in eachslice(_, dims=1)]
                reduce(hcat, _)
                exp.(unstandardise(_, log.(nbirds)))
            end
            # Wrap everything into a DataFrame for later summary in globalestimates.jl
            DataFrame(
                [atoll.unknown.str[pass], region.unknown.str[pass], species.unknown.str[pass], vec(μs), qs[1, :], qs[2, :]],
                [:atoll, :region, :species, :nbirds, :lower, :upper]
            )
        end
        Dict(Pair.(string.(thresholds), preds))
    end

# Save results for each threshold to an individual CSV
foreach(k -> CSV.write("$ROOT/results/data/countpreds_$k.csv", preds_target[k]), keys(preds_target))

# --- LOO CV --- #

# It seems that specifying the model as MvNormal confuses psis_loo ("1 data point")
# Respecify the likelihood as vectorized Normals
# y .~ Normal.(μ, σ)

# Regarding the "let's just put the model back into this script" ... 
# Not a big issue at all, but spelling out loospecial again seems nasty

loomodel = loospecial(inputs..., zlogn);

cv_res = psis_loo(loomodel, posterior.chains)
