# This script is part of the project associated with
# Article: Atolls are globally significant hubs for tropical seabirds
# Authors: Steibl S, Steiger S, Wegmann AS, Holmes ND, Young, HS, Carr P, Russell JC 
# Last edited: 2024-03-21

# --- MODEL MODULE --- #

"This module fits the count model, runs basic model diagnostics, summarises and exports the results."
module CountModel

export nothing

# --- WORKSPACE SETUP --- #

# Cmd args
# Need to check for "true" again because the shell script converts args to strings
load = Main.ARGS[1] == "false"
run_loocv = Main.ARGS[2] == "true"
run_sensitivity = Main.ARGS[3] == "true"

# Probabilistic programming
using Turing, TuringBenchmarking, ReverseDiff, ParetoSmooth, PosteriorStats
# Model optimization
using LazyArrays
# Statistics
using StatsFuns
# Manipulating input vectors
import StatsBase: denserank
# Covariance matrices
import LinearAlgebra: I
# Storing and working with data
using Chain, DataFrames, OrderedCollections
# Plotting
using StatsPlots
# Makie for ridgeline
import CairoMakie as CM
# Saving results
using Serialization, CSV
# Random seeds
using Random

# Load custom modules
include(joinpath(Main.ROOT, "julia", "src", "globalvars.jl"))
using .GlobalVariables
include(joinpath(Main.ROOT, "julia", "src", "countvars.jl"))
using .CountVariables
include(joinpath(Main.ROOT, "julia", "src", "utilities.jl"))
using .CustomUtilityFuns

# --- MODEL SETUP --- #

# Turing model
@model function model(
    r, s, n, PC, # input values
    idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes and array lengths
    y;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r), # generated indices and array
    pr=(μ_sxr=Ms, σ_sxr=[3, √2], μ_pxn=[0, 1], σ_pxn=[3, √2], σ2=[3, 2]) # prior, default prior is fallback
)

    # Priors for species × region
    μ_sxr ~ MvNormal(pr.μ_sxr, 1 * I)
    σ_sxr ~ filldist(InverseGamma(pr.σ_sxr...), Ns)
    z_sxr ~ filldist(Normal(), Ns, Nr)
    α_sxr = μ_sxr .+ σ_sxr .* z_sxr

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(pr.μ_pxn...), Nn, NPC)
    σ_pxn ~ filldist(InverseGamma(pr.σ_pxn...), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = μ_pxn[u_n, :] .+ σ_pxn[u_n, :] .* z_pxn[u_sn, :]

    # Prior for random error
    σ2 ~ InverseGamma(pr.σ2...)

    # Likelihood
    μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
    Σ = σ2 * I
    y ~ MvNormal(μ, Σ)

    # Generated quantities
    return (; y, α_sxr, β_pxn, σ2)
end

# Function to simulate samples from inputs
function simulate(params, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb; idx_sr=idx(s, r))
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
)


# Dictionary holding alternative prior settings for sensitivity analysis
dict_pr = Dict(
    "default" => (μ_sxr=Ms, σ_sxr=[3, √2], μ_pxn=[0, 1], σ_pxn=[3, √2], σ2=[3, 2]),
    # Differently informed prior, using mean instead of median
    "mean" => (μ_sxr=μs, σ_sxr=[3, √2], μ_pxn=[0, 1], σ_pxn=[3, √2], σ2=[3, 2]),
    # Using a much more weakly informed global prior in μ_sxr instead of a species-specific one
    "global" => (μ_sxr=fill(0, 37), σ_sxr=[3, √2], μ_pxn=[0, 1], σ_pxn=[3, √2], σ2=[3, 2]),
    # Narrow priors on variance parameters and μ_pxn, leaving μ_sxr at median
    "narrow" => (μ_sxr=Ms, σ_sxr=[3, √2 / 3], μ_pxn=[0, 1 / 3], σ_pxn=[3, √2 / 3], σ2=[3, 2 / 3]),
    # Wide priors on variance parameters and μ_pxn, leaving μ_sxr at median
    "wide" => (μ_sxr=Ms, σ_sxr=[3, √2 * 3], μ_pxn=[0, 1 * 3], σ_pxn=[3, √2 * 3], σ2=[3, 2 * 3]),
)

priorsettings = run_sensitivity ? collect(keys(dict_pr)) : ["default"]

# Define container for target predictions which will be updated inside the loop's local scope
preds_target = Dict{String,DataFrame}()

# Set autodiff to ReverseDiff{true}
Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

# Configure sampling
sampler = NUTS(1000, 0.95; max_depth=10)
nsamples = 10_000
nchains = 4
config = (sampler, MCMCThreads(), nsamples, nchains)

for priorsetting in priorsettings

    m = model(values(odict_inputs)..., zlogn; pr=dict_pr[priorsetting])

    # --- PRIOR PREDICTIVE CHECKS --- #

    @info "Count model: Prior predictive check for $priorsetting priors"

    # Set seed
    Random.seed!(42)

    # ModelSummary object holds chains, parameter names, and parameter samples
    prior = @chain begin
        sample(m, Prior(), 1000)
        ModelSummary(m, _)
    end

    let
        priorpreds = reduce(vcat, simulate(prior.samples, values(odict_inputs)...))
        histogram(zlogn, normalize=true, c=1)
        density!(priorpreds, normalize=true, c=:white, lw=3)
        density!(priorpreds, normalize=true, c=2, fillrange=0, fillalpha=0.2, legend=false)
        xlims!(-4, 4)
        title!("Prior predictive check, $priorsetting prior")
        xlabel!("z-standardised log-count")
        ylabel!("Density")
    end

    foreach(ext -> savefig(joinpath(Main.ROOT, "results", ext, "count", "prior_$(priorsetting)_$(Main.SUFFIX).$ext")), ["svg", "png"])

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

    # --- SAMPLING --- #

    # Notify user about sampling configuration
    @info """Count model: Sampling config for $priorsetting priors
    Sampler: $(string(sampler))
    Samples: $(nsamples)
    Chains: $(nchains)"""

    # Set seed
    Random.seed!(42)

    # ModelSummary object holds chains, parameter names, and parameter samples
    posterior = @chain begin
        if load
            try
                @info "Loading chains for $priorsetting prior."
                deserialize(joinpath(Main.ROOT, "results", "chains", "count_$(priorsetting).jls"))
            catch error
                @warn "Could not load chains, sampling from posterior instead."
                sample(m, config...)
            end
        else
            @info "Sampling from posterior."
            sample(m, config...)
        end
        ModelSummary(m, _)
    end

    diagnose(posterior.chains)

    !load && try
        path = joinpath(Main.ROOT, "results", "chains", "count_$(priorsetting)_$(Main.SUFFIX).jls")
        serialize(path, posterior.chains)
        @info "Count model: Chains saved to `$path`."
    catch error
        @warn "Count model: Writing chains failed with an $error."
    end

    # --- POSTERIOR PREDICTIVE CHECK --- #

    @info "Count model: Posterior predictive check for $priorsetting priors"

    preds_train = simulate(posterior.samples, values(odict_inputs)...)

    # Check if posterior similar to overall sample distribution
    # This check is of lower resolution than necessary for the research question
    # See species specific plots below
    let preds = reduce(vcat, preds_train)
        histogram(zlogn, normalize=true, c=1)
        density!(preds, normalize=true, c=:white, lw=3)
        density!(preds, normalize=true, c=2, fillrange=0, fillalpha=0.2, legend=false)
        xlims!(-4, 4)
    end

    # Overview plot of posterior predictions per species and atoll
    postpred_plot_view = let preds = exp.(unstandardise(reduce(hcat, preds_train), log.(nbirds)))
        sorted_species = sort(unique(species.known.str))
        # R reorders other vectors by region
        R = sortperm(region.known.num)

        # Enumerate sorted species to create one plot per species
        plots = map(enumerate(sorted_species)) do (idx, sp)
            # Create BitVector S for species subsetting
            S = species.known.num[R] .== idx
            obs = nbirds[R][S]
            # Calculate mean of predicted values for species S
            pred = vec(median(preds[R, :][S, :], dims=2))
            # Assemble plot for species S
            scatter(obs, markersize=2, msc=1, alpha=0.8, label="O", xformatter=_ -> "")
            scatter!(pred, markersize=2, msc=2, alpha=0.8, label="P")
            title!(replace(sp, "_" => " "), titlefontsize=12)
        end
        # Collect all plots in a dictionary (plotting all at once is a bit busy, maybe by nesting type or genus?)
        plot(plots..., layout=(8, 5), titlefontsize=8, size=(800, 1000))
    end

    foreach(ext -> savefig(joinpath(Main.ROOT, "results", ext, "count", "posterior_$(priorsetting)_$(Main.SUFFIX).$ext")), ["svg", "png"])

    # Create dictionary of species-wise posterior prediction plot, also showing uncertainty of estimates
    # Showing error bars makes the plot too busy for a single global view
    # Instead it is saved as a dictionary with one plot per species
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
                # Calculate mean of predicted values for species S
                pred = vec(median(preds[R, :][S, :], dims=2))
                # Calculate quantile of predicted values for species S
                lims = @chain preds[R, :][S, :] begin
                    map(x -> hdi(x, prob=0.9), eachslice(_, dims=1))
                    map(i -> abs.(getproperty.(Ref(_[i]), [:lower, :upper]) .- pred[i]), eachindex(_))
                    reduce(hcat, _)
                end
                # Assemble plot for species S
                scatter(atoll.known.str[R][S], obs, markersize=2.5, msc=2, label="O", permute=(:y, :x))
                scatter!(atoll.known.str[R][S], pred, yerror=(lims[1, :], lims[2, :]), markersize=2.5, msc=1, label="P", permute=(:y, :x))
                title!(replace(sp, "_" => " "), titlefontsize=12)
            end
            # Collect all plots in a dictionary (plotting all at once is a bit busy, maybe by nesting type or genus?)
            Dict(Pair.(sorted_species, plots))
        end

    # --- VALIDATION --- #

    # The data below consists of islands where only population ranges are known in the literature
    # We make predictions for these data and check how much of the posterior lies within these population ranges

    @info "Count model: Out-of-sample predictions for $priorsetting priors"

    preds_validation = let odict_oos_inputs = copy(odict_inputs)
        # Set names and values of inputs to be replaced in input dictionary
        names = ["region", "species", "PC", "species_within_nesting"]
        vals = [region.validation.num, species.validation.num, PC.validation, species_in_nesting.validation]
        setindex!.(Ref(odict_oos_inputs), vals, names)

        # Make predictions for validation data, transform to natural scale, and plot
        preds = @chain posterior.samples begin
            simulate(_, values(odict_oos_inputs)...)
            reduce(hcat, _)
            exp.(unstandardise(_, log.(nbirds)))
        end

        f = CM.Figure(size=(850, 1200))

        for (idx, sp) in enumerate(sort(unique(species.validation.str)))
            S = denserank(species.validation.num) .== idx
            x = minimum([3, size(preds[S, :], 1)])
            sppred = preds[S, :][1:x, :] # only species S, only x atolls
            lims = oos_lims[:, S][:, 1:x]

            q = maximum([quantile(slice, 0.8) for slice in eachslice(sppred, dims=1)])

            ax = CM.Axis(f[midx(idx, 4)...], title=replace(sp, "_" => " "), xlabel="Abundance (N)")

            for (i, v) in enumerate(eachslice(sppred, dims=1))
                CM.density!(v[v.<q],
                    strokecolor=:black,
                    strokewidth=0.5,
                    strokearound=true,
                    offset=-(i / q * 10),
                    color=:x,
                    colormap=blue_in_blue,
                    colorrange=(lims[1, i], lims[2, i])
                )
            end
            # Y axis represents three atolls, but labelling them would take up too much space on the view
            CM.hideydecorations!(ax)
            CM.hidespines!(ax, :t, :r, :l)
        end
        # Save plots to svg and png
        foreach(ext -> CM.save(joinpath(Main.ROOT, "results", ext, "count", "validation_$(priorsetting)_$(Main.SUFFIX).$ext"), f), ["svg", "png"])
    end

    # --- TARGET PREDICTIONS --- #

    # The data below are those atolls where no information on count data are present
    # We predict data for these atolls, using three different thresholds, and summarise it for export

    @info "Count model: Target predictions for $priorsetting model"

    global preds_target =
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
                    [getproperty.(Ref(hdi(slice; prob=0.95)), [:lower, :upper]) for slice in eachslice(_, dims=1)]
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
    foreach(k -> CSV.write(
            joinpath(Main.ROOT, "results", "data", "countpreds_$(k)_$(priorsetting)_$(Main.SUFFIX).csv"), select(preds_target[k], Not(:raw))
        ), keys(preds_target) # keys are k
    )

    # --- PSIS-LOO CV --- #

    # Specifying the model as MvNormal confuses psis_loo, the summary returns "sample size 1 data point"
    # Respecify the likelihood as vectorized Normals
    # y .~ Normal.(μ, σ)

    # Regarding the "let's just put the model back into this script" ... 
    # Not a big issue at all, but spelling out loospecial again seems nasty

    @model function broadcastmodel(
        r, s, n, PC, # input values
        idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes and array lengths
        y;
        Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r), # generated indices and array
        pr=(μ_sxr=Ms, σ_sxr=[3, √2], μ_pxn=[0, 1], σ_pxn=[3, √2], σ2=[3, 2])
    )

        # Priors for species × region
        μ_sxr ~ MvNormal(pr.μ_sxr, 1 * I)
        σ_sxr ~ filldist(InverseGamma(pr.σ_sxr...), Ns)
        z_sxr ~ filldist(Normal(), Ns, Nr)
        α_sxr = μ_sxr .+ σ_sxr .* z_sxr

        # Priors for nesting types × PCs
        μ_pxn ~ filldist(Normal(pr.μ_pxn...), Nn, NPC)
        σ_pxn ~ filldist(InverseGamma(pr.σ_pxn...), Nn, NPC)
        z_pxb ~ filldist(Normal(), Nb, NPC)
        z_pxg ~ filldist(Normal(), Ng, NPC)
        z_pxv ~ filldist(Normal(), Nv, NPC)
        z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
        β_pxn = μ_pxn[u_n, :] .+ σ_pxn[u_n, :] .* z_pxn[u_sn, :]

        # Prior for random error
        σ2 ~ InverseGamma(pr.σ2...)

        # Likelihood
        μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
        σ = sqrt(σ2)
        y .~ Normal.(μ, σ)

        # Generated quantities
        return (; y, α_sxr, β_pxn, σ2)
    end

    loomodel = broadcastmodel(values(odict_inputs)..., zlogn; pr=dict_pr[priorsetting])

    if run_loocv
        @info "Count model: Crossvalidation for $priorsetting priors"
        cv_res = psis_loo(loomodel, posterior.chains)
    else
        @warn "Count model: Skipping crossvalidation for $priorsetting priors"
    end

end

# --- GLOBAL ESTIMATES --- #

known = @chain pop_known begin
    select(_, [:atoll, :region, :species, :nbirds])
    transform(_, :nbirds => ByRow(x -> (lower=x, upper=x)) => AsTable)
end

glob = CSV.read(joinpath(Main.ROOT, "data", "birdlife_hbw_globalestimates_$(Main.SUFFIX).csv"), DataFrame)

DataFrames.transform!(pop_unknown, [:atoll, :species, :region] .=> denserank => x -> string("num_", x))

# Load data sets
p75, p80, p85 = CSV.read.(joinpath.(Main.ROOT, "results", "data", "countpreds_0." .* ["75", "8", "85"] .* "_default_$(Main.SUFFIX).csv"), DataFrame)
select!.([p75, p80, p85], Ref([:atoll, :region, :species, :median, :lower, :upper]))
rename!.([p75, p80, p85], :median => :nbirds)

full75, full80, full85 = [vcat(known, df) for df in [p75, p80, p85]]

tot75, tot80, tot85 = popsum.([full75, full80, full85])
leftjoin!.([tot75, tot80, tot85], Ref(glob), on=:species)

for n in [:nbirds, :lower, :upper], df in [tot75, tot80, tot85]
    transform!(df, [:species, n, :birdlife_min, :birdlife_max, :HBW, :Otero] => ByRow(calcratio) => string("ratio_", n))
end
select!.([tot75, tot80, tot85], Ref(Not(:birdlife_min, :birdlife_max, :HBW, :Otero)))

histogram((tot75.ratio_nbirds .- tot85.ratio_nbirds) .* 100, title="Pop on atoll ratio diff .75 - .85", label=:none)
xticks!(-1:3, string.(-1:3) .* "%")

CSV.write(joinpath(Main.ROOT, "results", "data", "pred_and_obs_atolls_$(Main.SUFFIX).csv"), full80)

# Summarise raw samples before further summaries

full_with_raw = @chain known begin
    transform(_, :nbirds => ByRow(x -> fill(x, Int(nsamples*nchains))) => :raw)
    rename(_, :nbirds => :median)
    select(_, :atoll, :region, :species, :raw, :median, :lower, :upper)
    vcat(preds_target["0.8"], _)
end

summary_specieswise = let df = full_with_raw
    out = DataFrame()
    for sp in unique(df.species)
        tmp = @chain df begin
            subset(_, :species => ByRow(x -> x == sp))
            getproperty(_, :raw)
            sum(_)
            Dict("species" => sp, "median" => median(_), "lower" => quantile(_, 0.025), "upper" => quantile(_, 0.975))
        end
        append!(out, DataFrame(tmp))
    end
    out
end

summary_atollwise = let df = full_with_raw
    out = DataFrame()
    for atoll in unique(df.atoll)
        tmp = @chain df begin
            subset(_, :atoll => ByRow(x -> x == atoll))
            getproperty(_, :raw)
            sum(_)
            Dict("atoll" => atoll, "median" => median(_), "lower" => quantile(_, 0.025), "upper" => quantile(_, 0.975))
        end
        append!(out, DataFrame(tmp))
    end
    out
end

CSV.write(joinpath(Main.ROOT, "results", "data", "summary_specieswise_$(Main.SUFFIX).csv"), summary_specieswise)
CSV.write(joinpath(Main.ROOT, "results", "data", "summary_atollwise_$(Main.SUFFIX).csv"), summary_atollwise)

# --- SENSITIVITY ANALYSIS --- # 

# for 0.8 cutoffs

dict_sensitivity = @chain begin
    map(priorsettings) do priorsetting
        Pair(priorsetting, CSV.read(joinpath(Main.ROOT, "results", "data", "countpreds_0.8_$(priorsetting)_$(Main.SUFFIX).csv"), DataFrame))
    end
    Dict(_)
end

select!.(values(dict_sensitivity), Ref([:atoll, :region, :species, :median, :lower, :upper]))
rename!.(values(dict_sensitivity), :median => :nbirds)
[dict_sensitivity[k] = vcat(known, v) for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))]
[dict_sensitivity[k] = popsum(v) for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))]
leftjoin!.(values(dict_sensitivity), Ref(glob), on=:species)

for n in [:nbirds, :lower, :upper], df in values(dict_sensitivity)
    transform!(df, [:species, n, :birdlife_min, :birdlife_max, :HBW, :Otero] => ByRow(calcratio) => string("ratio_", n))
end

for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))
    dict_sensitivity[k] = insertcols(v, 1, :prior => k)
end

# Compare global predictions for different prior settings
df_sens = reduce(vcat, values(dict_sensitivity))

scatter(df_sens.ratio_nbirds, df_sens.species, group=df_sens.prior, alpha=0.75, markershape=:x)
xlabel!("Global population on atolls")
xticks!(0:0.5:2, string.(Int64.(collect(0:0.5:2) .* 100), "%"))
yticks!(eachindex(unique(df_sens.species)), unique(df_sens.species), size=(600, 600), tickfontsize=7)

foreach(ext -> savefig(joinpath(Main.ROOT, "results", ext, "count", "sensitivity_count.$ext")), ["svg", "png"])

end