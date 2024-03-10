# This script is part of the project associated with
# Article: Atolls are globally significant hubs for tropical seabirds
# Authors: Steibl S, Steiger S, Wegmann AS, Holmes ND, Young, HS, Carr P, Russell JC 
# Last edited: 2024-03-10

# --- MODEL MODULE --- #

"This module fits the presence model, runs basic model diagnostics, summarises and exports the results."
module PresenceModel

export nothing

# --- WORKSPACE SETUP --- #

# Cmd args
# Need to check for "true" again because the shell script converts args to strings
load = Main.ARGS[1] == "false"
run_loocv = Main.ARGS[2] == "true"
run_sensitivity = Main.ARGS[3] == "true"

# Probabilistic programming
using Turing, ParetoSmooth, ReverseDiff
# Benchmarking
using TuringBenchmarking
# Model speed optimization
using LazyArrays
# Statistics
using StatsFuns
# Variance matrices
import LinearAlgebra: I
# Working with tabular data
using Chain, DataFrames
# Plotting
using StatsPlots
# Saving results and logging
using Serialization, CSV, Dates, OrderedCollections
# Random seeds
using Random

# Load custom modules
include(joinpath(Main.ROOT, "julia", "src", "presencevars.jl"))
using .PresenceVariables
include(joinpath(Main.ROOT, "julia", "src", "utilities.jl"))
using .CustomUtilityFuns

# --- MODEL SETUP --- #

# Turing model
@model function model(
    r, s, n, PC, # inputs
    idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes, lengths, etc
    y; # outcome
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r),
    pr=(α_sxr=[0, 1], μ_pxn=[0, 0.2], σ_pxn=[3, 0.5])
)

    # Priors for species × region
    α_sxr ~ filldist(Normal(pr.α_sxr...), Ns * Nr)

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(pr.μ_pxn...), Nn, NPC)
    σ_pxn ~ filldist(InverseGamma(pr.σ_pxn...), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = @. μ_pxn[u_n, :] + σ_pxn[u_n, :] * z_pxn[u_sn, :]

    # Likelihood
    v = α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2)
    y .~ BernoulliLogit.(v)

    # Generated quantities
    return (; α_sxr, β_pxn)
end

# Function to simulate samples from inputs to model
function simulate(params, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb; idx_sr=idx(s, r))
    map(params) do p
        β = sum(p.β_pxn[idx_sn, :] .* X, dims=2)
        @. logistic(p.α_sxr[idx_sr] + β)
    end
end

# Dictionary holding inputs
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

# Dictionary holding alternative prior settings for sensitivity analysis
dict_pr = Dict(
    "default" => (α_sxr=[0, 1], μ_pxn=[0, 0.2], σ_pxn=[3, 0.5]),
    # Differently informed prior, using mean instead of median
    "narrow" => (α_sxr=[0, 1 / 3], μ_pxn=[0, 0.2 / 3], σ_pxn=[3, 0.5 / 3]),
    # Wide priors on variance parameters and μ_pxn, leaving μ_sxr at median
    "wide" => (α_sxr=[0, 1 * 3], μ_pxn=[0, 0.2 * 3], σ_pxn=[3, 0.5 * 3]),
)

priorsettings = run_sensitivity ? collect(keys(dict_pr)) : ["default"]

for priorsetting in priorsettings

    m = model(values(odict_inputs)..., presence; pr=dict_pr[priorsetting])

    # --- PRIOR PREDICTIVE CHECK --- #

    @info "Presence model: Prior predictive check for $priorsetting priors"

    # Set seed
    Random.seed!(42)

    # ModelSummary object holds chains, parameter names, and parameter samples
    prior = @chain begin
        sample(m, Prior(), 1000)
        ModelSummary(m, _)
    end

    let
        priorpreds = reduce(vcat, simulate(prior.samples, values(odict_inputs)...))
        density(priorpreds, normalize=true, c=2, fillrange=0, fillalpha=0.2, legend=false)
        title!("Prior predictive check, $priorsetting prior")
        ylabel!("Density"), xlabel!("Probability of presence")
    end

    foreach(ext -> savefig(joinpath(Main.ROOT, "results", ext, "presence", "prior_$priorsetting.$ext")), ["svg", "png"])

    # --- MODEL CONFIG --- #

    # We tested several autodiff methods to check which one is fastest for our model
    # ReverseDiff with caching enabled is the fastest choice.
    benchmark = false # do not rerun unless specifically requested

    if benchmark
        @info "Presence model: Benchmarking model for $priorsetting priors"

        backends = [
            Turing.Essential.ForwardDiffAD{0}(),
            Turing.Essential.ReverseDiffAD{false}(),
            Turing.Essential.ReverseDiffAD{true}()
        ]

        TuringBenchmarking.run(TuringBenchmarking.make_turing_suite(m, adbackends=backends);)
    end

    # Set autodiff to ReverseDiff{true}
    Turing.setadbackend(:reversediff)
    Turing.setrdcache(true)

    # Configure sampling
    sampler = NUTS(1000, 0.95; max_depth=10)
    nsamples = 10_000
    nchains = 4
    config = (sampler, MCMCThreads(), nsamples, nchains)

    @info """Presence model: Sampling config for $priorsetting priors
    Sampler: $(string(sampler))
    Samples: $(nsamples)
    Chains: $(nchains)"""

    # Set seed
    Random.seed!(42)

    posterior = @chain begin
        if load
            try
                @info "Loading chains for $priorsetting priors."
                deserialize(joinpath(Main.ROOT, "results", "chains", "presence_$(priorsetting).jls"))
            catch error
                @warn "Loading failed with an $error, sampling from posterior instead."
                sample(m, config...) # Sample from model
            end
        else
            @info "Sampling from posterior."
            sample(m, config...) # Sample from model
        end
        ModelSummary(m, _)
    end

    diagnose(posterior.chains)


    !load && try
    path = joinpath(Main.ROOT, "results", "chains", "presence_$(priorsetting)_$(Main.SUFFIX).jls")
        serialize(path, posterior.chains)
        @info "Presence model: Chains saved to `$path`."
    catch error
        @warn "Presence model: Writing chains failed with an $error."
    end

    # --- POSTERIOR PREDICTIVE CHECK --- #

    @info "Presence model: Posterior predictive check for $priorsetting priors"

    postpcplots = let preds = reduce(hcat, simulate(posterior.samples, values(odict_inputs)...))
        # Iterate over species and create plots with posterior predictions
        map(enumerate(unique(species.known.num))) do (idx, sp)
            # Predicted values for species sp
            pred_x = vec(preds[species.known.num.==sp, :])
            # Observed values for species sp
            obs_x = presence[species.known.num.==sp]
            # Assemble plot for species sp
            histogram(obs_x, normalize=true, alpha=0.5, lc=:transparent, bins=10, label="O")
            density!(pred_x, fillrange=0, fillalpha=0.2, normalize=true, alpha=0.5, lc=:transparent, yticks=:none, label="P")
            xticks!(0:0.5:1, string.(0:0.5:1))
            title!(unique(species.known.str)[idx], titlefontsize=8)
            vline!([0.8], c=:black, ls=:dash, label=:none)
        end
    end

    plot(postpcplots..., titlefontsize=9, size=(800, 1200), layout=(8, 5))
    foreach(ext -> savefig("$(Main.ROOT)/results/$ext/presence/posterior_$(priorsetting)_$(Main.SUFFIX).$ext"), ["svg", "png"])

    # --- TARGET PREDICTIONS --- #

    @info "Presence model: Target predictions for $priorsetting model"

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
        plots = map(enumerate(unique(df.species))) do (idx, sp)
            I = df.species .== sp
            scatter(df.mean[I], df.atoll[I], c=ifelse.(df.mean[I] .> 0.8, :red, :black), ms=2, xerror=(abs.(df.lower005[I] .- df.mean[I]), abs.(df.upper095[I] .- df.mean[I])))
            title!(unique(species.unknown.str)[idx], titlefontsize=8)
            xlims!(0, 1)
            vline!([0.8], c=:black, ls=:dash, legend=:none)
        end
        Dict(Pair.(unique(species.unknown.str), plots))
    end

    # plot(dict_unknown_pred_plots..., titlefontsize=9, size=(800, 1200), layout=(8, 5))
    # TODO plot by nesting type or genus

    # Save predictions for unknown atolls to CSV
    try
        path = joinpath(Main.ROOT, "results", "data", "presencepreds_$(priorsetting)_$(Main.SUFFIX).csv")
        CSV.write(path, df_preds)
        @info "Presence model: Saved predictions to `$path`."
    catch error
        @warn "Presence model: Writing predictions failed with an $error."
    end

    # --- PSIS-LOO CV --- #
    if run_loocv
        @info "Presence model: Crossvalidation for $priorsetting priors"
        cv_res = psis_loo(m, posterior.chains)
        cv_res
        # - no overfit
        # - out of sample performance near in-sample performance (gmpd 0.76)
        # - not many outliers in the p_eff plot (the outliers are logical => Clipperton, Ant (PCs?))
        # - in line with posterior predictive check
    else
        @warn "Presence model: Skipping crossvalidation for $priorsetting priors"
    end

end

end

# --- SENSITIVITY MODULE --- #

"This module runs the sensitivity analyses for the presence model and exports the results."
module PresenceSensitivity

using DataFrames, Chain, CSV, StatsPlots

# Helper to unstack dataframe
my_unstack(df) = unstack(df[:, [:atoll, :species, :mean]], :species, :mean)

# Only running cutoff sensitivity checks for default priors here
preds_default = my_unstack(CSV.read(joinpath(Main.ROOT, "results", "data", "presencepreds_default_$(Main.SUFFIX).csv"), DataFrame))

plot(
    map([0.7:0.05:0.9;]) do threshold
        heatmap(Matrix(preds_default[!, 2:end] .> threshold), title="Threshold: $threshold")
    end...,
    layout=(3, 2),
    size=(600, 800),
    plot_title="Species included per cutoff",
)
ylabel!("Atoll")
xlabel!("Species")

foreach(ext -> savefig(joinpath(Main.ROOT, "results", ext, "presence", "sensitivity_$(Main.SUFFIX).$ext")), ["svg", "png"])

end