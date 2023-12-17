# --- WORKSPACE SETUP --- #

# Probabilistic programming
using Turing, TuringBenchmarking, ReverseDiff, ParetoSmooth
# Model speed optimization
using LazyArrays
# Statistics
using StatsFuns
import LinearAlgebra: I
# Working with tabular data
using Chain, DataFrames
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

# Set seed
Random.seed!(42)

# Benchmark model?
benchmark = false

# Run the sampler?
run = isempty(ARGS) ? false : ARGS[1] == "true"

model = mMc
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
chainpath = "count_default.jls"

standardise(x) = (x .- mean(x)) ./ std(x)

zlogn = standardise(log.(nbirds))

# --- PRIOR PREDICTIVE CHECKS --- #

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
benchmark && let
    backends = [
        Turing.Essential.ForwardDiffAD{0}(),
        Turing.Essential.ReverseDiffAD{false}(),
        Turing.Essential.ReverseDiffAD{true}()
    ]

    TuringBenchmarking.run(TuringBenchmarking.make_turing_suite(model(inputs..., presence), adbackends=backends);)
end

# Sample from model unless a saved chain should be used
if !run
    chain = deserialize("$ROOT/results/chains/$chainpath")
else
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
    Chains: $(nchains)
    """

    # Discarding samples is unnecessary after NUTS tuning
    @info "🚀 Starting sampling: $(Dates.now())"
    posterior = @chain begin
        sample(m, config...)
        ModelSummary(m, _)
    end

    serialize("$ROOT/results/chains/$chainpath", chain)
    @info "💾 Chain saved to '$ROOT/results/chains/$chainpath'."
end;

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
    # Create matrix of posterior predictions, rows X is species X
    preds = reduce(hcat, posteriorpreds)
    # Enumerate species keys to create one plot per species
    plots = map(enumerate(keys(odict_species))) do (index, species)
        # Get observed values for species X
        obs = zlogn[num_species_known.==index]
        # Calculate mean of predicted values for species X
        pred = vec(mean(preds[num_species_known.==index, :], dims=2))
        # Calculate quantile of predicted values for species X
        qlims = @chain preds[num_species_known.==index, :] begin
            map(x -> quantile(x, [0.05, 0.95]), eachslice(_, dims=1))
            map(i -> preds[i] .- _[i], eachindex(_))
            transpose(reduce(hcat, _))
        end
        # Assemble plot for species X
        scatter(pred, yerror=qlims, markersize=2.5, msc=1, label="P")
        scatter!(obs, markersize=2.5, msc=2, label="O", xticks=:none)
        title!(species, titlefontsize=12)
    end
    # Collect all plots in a dictionary (plotting all at once is a bit busy, maybe by nesting type or genus?)
    Dict(Pair.(keys(odict_species), plots))
end

# TODO add atoll names as x label, then flip xy

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

# TODO move to utils
between(x, lower, upper) = lower ≤ x && x ≤ upper

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

loomodel = loospecial(inputs..., zlogn);

cv_res = psis_loo(loomodel, chain)
