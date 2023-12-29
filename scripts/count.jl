module CountModel

export nothing

# --- WORKSPACE SETUP --- #

# Probabilistic programming
using Turing, TuringBenchmarking, ReverseDiff, ParetoSmooth, PosteriorStats
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
using Serialization, CSV, Dates
# Random seeds
using Random

# Path
const ROOT = dirname(Base.active_project())

# Load custom modules
include("$ROOT/src/count.jl")
using .CountVariables
include("$ROOT/src/utilities.jl")
using .CustomUtilityFuns

# --- MODEL SETUP --- #

# Turing model
@model function model(
    r, s, n, PC, # input values
    idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes
    ssp,
    y;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r) # generated indices
)

    # Priors for species × region
    μ_sxr ~ MvNormal(ssp, 1 * I)
    σ_sxr ~ filldist(InverseGamma(3, 2), Ns)
    z_sxr ~ filldist(Normal(), Ns, Nr)
    α_sxr = μ_sxr .+ σ_sxr .* z_sxr

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, 1), Nn, NPC)
    σ_pxn ~ filldist(InverseGamma(3, 2), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = μ_pxn[u_n, :] .+ σ_pxn[u_n, :] .* z_pxn[u_sn, :]

    # Prior for random error
    σ2 ~ InverseGamma(3, 2)

    # Likelihood
    μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
    Σ = σ2 * I
    y ~ MvNormal(μ, Σ)

    # Generated quantities
    return (; y, α_sxr, β_pxn, σ2)
end

# Function to simulate samples from inputs
function simulate(params, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb, ssp; idx_sr=idx(s, r))
    map(enumerate(params)) do (i, param)
        μ = param.α_sxr[idx_sr] + sum(param.β_pxn[idx_sn, :] .* X, dims=2)
        Random.seed!(i)
        rand.(Normal.(μ, param.σ2))
    end
end

# Set inputs
odict_inputs = OrderedDict(
    "region" => region.known.num,
    "species" => species.known.num,
    "nesting" => nesting.known.num,
    "PC" => PC.known,
    "species_within_nesting" => species_in_nesting.known,
    "unique_nesting" => nesting.levels,
    "unique_species_within_nesting" => species_in_nesting.levels,
    "n_burow" => nspecies.burrow,
    "n_ground" => nspecies.ground,
    "n_vegetation" => nspecies.vegetation,
    "ssp" => median_species_zlogn,
)

# If not loading a chain, save results to path below
chainpath = "count.jls"

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
    xlims!(-10, 10)
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

# Set autodiff to ReverseDiff{true}
Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

# Configure sampling
sampler = NUTS(1000, 0.95; max_depth=10)
nsamples = 10_000
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
    path = "$ROOT/results/chains/$chainpath"
    serialize(path, posterior.chains)
    @info "💾 Chain saved to `$path`."
catch error
    @warn "Writing chains failed with an $error."
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
        # R reorders other vectors by region
        R = sortperm(region.known.num)

        # Enumerate sorted species to create one plot per species
        plots = map(enumerate(sorted_species)) do (idx, sp)
            # Create BitVector S for species subsetting
            S = species.known.num[R] .== idx
            obs = nbirds[R][S]
            # Calculate mean of predicted values for species X
            pred = vec(median(preds[R, :][S, :], dims=2))
            # Calculate quantile of predicted values for species X
            qlims = @chain preds[R, :][S, :] begin
                map(x -> hdi(x, prob=0.9), eachslice(_, dims=1))
                map(i -> abs.(getproperty.(Ref(_[i]), [:lower, :upper]) .- pred[i]), eachindex(_))
                reduce(hcat, _)
            end
            # Assemble plot for species X
            scatter(atoll.known.str[R][S], pred, yerror=(qlims[1, :], qlims[2, :]), markersize=2.5, msc=1, label="P", permute=(:y, :x))
            scatter!(atoll.known.str[R][S], obs, markersize=2.5, msc=2, label="O", permute=(:y, :x))
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
    mdn = @chain posterior.samples begin
        simulate(_, values(odict_oos_inputs)...)
        reduce(hcat, _)
        exp.(unstandardise(_, log.(nbirds)))
        vec(median(_, dims=2))
    end
    pass = vec(Base.isbetween.(oos_lims[1, :], mdn, oos_lims[2, :]))
    DataFrame([species.validation.str, oos_lims[1, :], mdn, oos_lims[2, :], pass], [:species, :lower, :pred, :upper, :pass])
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
            raw = [exp.(unstandardise(slice, log.(nbirds))) for slice in eachslice(samples, dims=1)]
            mdn = exp.(unstandardise(median(samples, dims=2), log.(nbirds)))
            qs = @chain samples begin
                [getproperty.(Ref(hdi(slice; prob=0.75)), [:lower, :upper]) for slice in eachslice(_, dims=1)]
                reduce(hcat, _)
                exp.(unstandardise(_, log.(nbirds)))
            end
            # Wrap everything into a DataFrame for later summary in globalestimates.jl
            DataFrame(
                [atoll.unknown.str[pass], region.unknown.str[pass], species.unknown.str[pass], raw, vec(mdn), qs[1, :], qs[2, :]],
                [:atoll, :region, :species, :raw, :median, :lower, :upper]
            )
        end
        Dict(Pair.(string.(thresholds), preds))
    end

plots_preds_target = map(eachrow(preds_target["0.8"])) do row
    density(row.raw[row.raw.<quantile(row.raw, 0.95)], title="$(row.species) on $(row.atoll) ($(row.region))", label=:none, fillrange=0, fillalpha=0.2)
    vline!([median(row.raw)], c=:red, lw=2, label="med=$(round(median(row.raw), digits=0))", alpha=0.5)
    vline!([row.lower, row.upper], c=:black, ls=:dash, alpha=0.5, label="$(round(row.lower, digits=0))-$(round(row.upper, digits=0))")
end

# Save results for each threshold to an individual CSV
foreach(k -> CSV.write("$ROOT/results/data/countpreds_$k.csv", select(preds_target[k], Not(:raw))), keys(preds_target))

# --- LOO CV --- #

# It seems that specifying the model as MvNormal confuses psis_loo ("1 data point")
# Respecify the likelihood as vectorized Normals
# y .~ Normal.(μ, σ)

# Regarding the "let's just put the model back into this script" ... 
# Not a big issue at all, but spelling out loospecial again seems nasty

@model function broadcastmodel(
    r, s, n, PC,
    idx_sn, u_n, u_sn, Nv, Ng, Nb,
    y;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    μ_sxr ~ filldist(Normal(0, 0.2), Ns)
    τ_sxr ~ filldist(InverseGamma(3, 0.2), Ns)
    z_sxr ~ filldist(Normal(), Ns, Nr)
    α_sxr = μ_sxr .+ τ_sxr .* z_sxr

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, 0.1), Nn, NPC)
    τ_pxn ~ filldist(InverseGamma(3, 0.1), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = μ_pxn[u_n, :] .+ τ_pxn[u_n, :] .* z_pxn[u_sn, :]

    # Prior for random error
    σ2 ~ InverseGamma(3, 0.2)

    # Likelihood
    μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
    σ = sqrt(σ2)
    y .~ Normal.(μ, σ)

    # Generated quantities
    return (; y, α_sxr, β_pxn, σ2)
end

loomodel = broadcastmodel(inputs..., zlogn)

cv_res = psis_loo(loomodel, posterior.chains)

end
