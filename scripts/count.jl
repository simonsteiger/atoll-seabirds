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

# Prior settings can be set from command line
# σₚ, θₚ = 1, 2

# if !run
#     @info "Loading chain, no model fit."
# elseif isempty(ARGS) || ARGS[2] == "default"
#     @info "Fitting model with default priors: σₚ=$σₚ, θₚ=$θₚ."
# elseif all(ARGS[2] .!= ["narrow", "wide"])
#     throw("Unknown prior setting: '$(ARGS[2])'. Pass nothing or one of 'default', 'narrow', 'wide'.")
# else
#     σₚ, θₚ = ARGS[2] == "wide" ? [σₚ, θₚ] .* 3 : [σₚ, θₚ] .* 1 / 3
#     @info "Fitting model with $(ARGS[2]) priors: σₚ=$(round(σₚ, digits=2)), θₚ=$(round(θₚ, digits=2))."
# end

if ARGS[2] == "S" # simplest model
    # model = m1p
    # inputs = (; num_region, num_species)
elseif isempty(ARGS) || ARGS[2] == "M" # most complex model
    model = mMc
    # What does the auto named tuple notation (; ...) do with splatted vectors?
    inputs = (num_region_known, num_species_known, num_nesting_known, PC_known, num_species_within_nesting_known, unique_nesting_known, unique_species_within_nesting_known, count_species_by_nesting...,)
elseif isempty(ARGS) || ARGS[2] == "L" # large model
    # model = mLp
    # What does the auto named tuple notation (; ...) do with splatted vectors?
    # inputs = (num_region, num_species, num_nesting, PC, num_species_within_nesting, unique_nesting, unique_species_within_nesting, count_species_by_nesting...,)
elseif ARGS[2] == "XL"
    # model = mXLp
    # What does the auto named tuple notation (; ...) do with splatted vectors?
    # inputs = (num_atoll, num_region, num_species, num_nesting, PC, num_species_within_nesting, unique_nesting, unique_species_within_nesting, count_species_by_nesting...,)
end

PRIORSUFFIX = isempty(ARGS) ? "default" : ARGS[2]

# If not loading a chain, save results to path below
chainpath = "count.jls"

standardise(x) = (x .- mean(x)) ./ std(x)

zlogn = standardise(log.(nbirds))

# --- PRIOR PREDICTIVE CHECKS --- #

chain_prior = sample(model(inputs..., zlogn), Prior(), 1000);

θ_prior = let
    peaceful_generated_quantities(model(inputs..., zlogn), chain_prior);
    k = keys(θ_prior[1])
    vec(getsamples(θ_prior, k...))
end

let 
    priorsamples = reduce(vcat, simulate(θ_prior, inputs...))
    histogram(zlogn, normalize=true, c=1)
    density!(priorsamples, normalize=true, c=:white, lw=3)
    density!(priorsamples, normalize=true, c=2, fillrange=0, lw=1.5, fillalpha=0.2, legend=false)
end 
# not sure if necessary to show both observed and prior for prior predictive check
# I guess it's not bad on a second thought though – we don't want probability mass in useless places

# Benchmark different backends to find out which is fastest
let adbackends = [:forwarddiff, :reversediff, :reversediff_compiled]
    benchmark && benchmark_model(model; check=true, adbackends=adbackends)
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
    ndiscard = 200

    @info """Sampler: $(string(sampler))
    Samples: $(nsamples)
    Chains: $(nchains)
    Discard: $(ndiscard)
    """

    @info "🚀 Starting sampling: $(Dates.now())"
    chain = sample(model(inputs..., zlogn), sampler, MCMCThreads(), nsamples, nchains; discard_initial=ndiscard)

    serialize("$ROOT/results/chains/$chainpath", chain)
    @info "💾 Chain saved to '$ROOT/results/chains/$chainpath'."
end;

# --- POSTERIOR PREDICTIVE CHECK --- #

θ_post = let
    gqs = peaceful_generated_quantities(model(inputs..., zlogn), chain)
    k = keys(θ_post[1])
    vec(getsamples(gqs, k...))
end

let 
    postsamples = reduce(hcat, simulate(θ_post, inputs...))

    plots = map(enumerate(unique(num_species_known))) do (index, species)
        pred_x = vec(mean(postsamples[num_species_known.==species, :], dims=2))
        obs_x = zlogn[num_species_known.==species]
        scatter(obs_x, markersize=2.5, msc=1, label="O")
        scatter!(pred_x, markersize=2.5, msc=2, yticks=:none, label="P")
        #title!(unique(str_species)[index], titlefontsize=8)
    end

    plot(plots..., titlefontsize=9, size=(800,1200), layout=(8,5))
end

# Goal here is scatter plots that also show uncertainty in estimates – maybe even violins? Too busy...

predictions = 
    let thresholds = 0.75:0.05:0.85 
        @chain thresholds begin
            map(_) do threshold
                pass = ppres .> threshold
                unknown_inputs = (num_region_unknown[pass], num_species_unknown[pass], num_nesting_known, PC_unknown[pass, :], num_species_within_nesting_unknown[pass], unique_nesting_known, unique_species_within_nesting_known, count_species_by_nesting...,)
                samples = reduce(hcat, simulate(θ_post, unknown_inputs...))
            end
        Dict(Pair.(thresholds, _))
    end
end


countpreds_unknown = @chain begin
    predictcount(
        α,
        β,
        σ2,
        num_species_within_nesting_unknown[threshold],
        num_species_unknown[threshold],
        num_region_unknown[threshold],
        PC_unknown[threshold, :],
    )
    reduce(hcat, _)
    Matrix{Float64}(_)
end

between(x, lower, upper) = lower ≤ x && x ≤ upper

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



countpreds_known = @chain begin
    predictcount(
        α,
        β,
        σ2,
        num_species_within_nesting_known,
        num_species_known,
        num_region_known,
        PC_known,
    )
    reduce(hcat, _)
    Matrix{Float64}(_)
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

# Overall posterior predictive checks
hist_overall_postpc = histogram(log.(nbirds), normalize=true, lw=0.5, lc=:white, label="Observed")
histogram!(vec(countpreds_known), normalize=true, c=:white, lw=3, label=false)
histogram!(vec(countpreds_known), normalize=true, c=2, lw=1.5, label="Predicted")
title!("Posterior predictive check"); xlabel!("log(Count)"); ylabel!("Density")

# Gather prior and posterior predictive checks into single plot
grid_pppc = plot(hist_overall_priorpc, hist_overall_postpc, size=(800, 500))
xlims!(-30, 30)
display(grid_pppc)
savefig("$ROOT/results/svg/validation_count_pppc.svg")

# Species-specific checks
hist_species_postpc = map(enumerate(unique(num_species_known))) do (index, species)
    pred_x = avg_preds_known[num_species_known.==species]
    obs_x = log.(nbirds)[num_species_known.==species]
    scatter(eachindex(pred_x), obs_x, markersize=2.5, msc=1, label="O")
    scatter!(eachindex(pred_x), pred_x, markersize=2.5, msc=2, label="P", xformatter=_->"")
    title!(unique(str_species_known)[index])
end

grid_species_ppc = plot(hist_species_postpc..., layout=(8,5), size=(1000, 1400), titlefontsize=9)
savefig("$ROOT/results/svg/validation_count_postpc_species.svg")


# --- LOO CV --- #

# It seems that specifying the model as MvNormal breaks psis_loo ("1 data point")
# Respecify the likelihood as vectorized Normals
# y .~ Normal.(μ, σ)
# Can still use the fit chain though because formulas are identical.

@model function loospecial(
    r, s, n, PC, y,
    idx_sn, u_n, u_sn, Nv, Ng, Nb;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    μ_sxr ~ filldist(Normal(0, σₚ), Ns)
    τ_sxr ~ filldist(InverseGamma(3, θₚ), Ns)
    z_sxr ~ filldist(Normal(), Ns, Nr)
    α_sxr = μ_sxr .+ τ_sxr .* z_sxr

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, σₚ), Nn, NPC)
    τ_pxn ~ filldist(InverseGamma(3, θₚ), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = μ_pxn[u_n, :] .+ τ_pxn[u_n, :] .* z_pxn[u_sn, :]

    # Prior for random error
    σ2 ~ InverseGamma(3, θₚ)

    # Likelihood
    μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
    y .~ Normal.(μ, σ2)

    # Generated quantities
    return (; y, α_sxr, β_pxn)
end;

loomodel = loospecial(
    num_region_known,
    num_species_known,
    num_nesting_known,
    PC_known,
    log.(nbirds),
    num_species_within_nesting_known,
    unique_nesting_known,
    unique_species_within_nesting_known,
    count_species_by_nesting...,
);

cv_res = psis_loo(loomodel, chain)

# Variance priors: Exponential(1).^2  #
# ----------------------------------- #
# Warning: Some Pareto k values are extremely high (>1). PSIS will not produce consistent estimates.
# Results of PSIS-LOO-CV with 8000 Monte Carlo samples and 861 data points. Total Monte Carlo SE of 0.45.
# ┌───────────┬──────────┬──────────┬───────┬─────────┐
# │           │    total │ se_total │  mean │ se_mean │
# ├───────────┼──────────┼──────────┼───────┼─────────┤
# │   cv_elpd │ -1927.46 │    26.03 │ -2.24 │    0.03 │
# │ naive_lpd │ -1787.33 │    18.71 │ -2.08 │    0.02 │
# │     p_eff │   140.13 │     8.97 │  0.16 │    0.01 │
# └───────────┴──────────┴──────────┴───────┴─────────┘

# Variance priors: InverseGamma(3, 0.1)  #
# -------------------------------------- #
# No warning
# Results of PSIS-LOO-CV with 8000 Monte Carlo samples and 861 data points. Total Monte Carlo SE of 0.11.
# ┌───────────┬──────────┬──────────┬───────┬─────────┐
# │           │    total │ se_total │  mean │ se_mean │
# ├───────────┼──────────┼──────────┼───────┼─────────┤
# │   cv_elpd │ -2221.51 │     4.99 │ -2.58 │    0.01 │
# │ naive_lpd │ -2215.20 │     4.74 │ -2.57 │    0.01 │
# │     p_eff │     6.31 │     0.33 │  0.01 │    0.00 │
# └───────────┴──────────┴──────────┴───────┴─────────┘
