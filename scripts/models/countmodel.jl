# Prior settings can be set from command line
σₚ, λₚ = 1, 1

if isempty(ARGS) || ARGS[1] == "default"
    @info "Fitting model with default priors: σₚ=$σₚ, λₚ=$λₚ."
elseif all(ARGS[1] .!= ["narrow", "wide"])
    throw("Unknown prior setting: '$(ARGS[1])'. Pass nothing or one of 'default', 'narrow', 'wide'.")
else
    σₚ, λₚ = ARGS[1] == "wide" ? [σₚ, λₚ] .* 3 : [σₚ, λₚ] .* 1 / 3
    @info "Fitting model with $(ARGS[1]) priors: σₚ=$(round(σₚ, digits=2)), λₚ=$(round(λₚ, digits=2))."
end

PRIORSUFFIX = isempty(ARGS) ? "default" : ARGS[1]

const ROOT = dirname(Base.active_project())

# Probabilistic programming
using Turing, TuringBenchmarking, ReverseDiff
# Model speed optimization
using LazyArrays
# Statistics
using StatsFuns, LinearAlgebra
# Working with tabular data
using Chain, DataFrames
# Plotting
using StatsPlots
# Saving results and logging
using Serialization, CSV, Dates, Markdown
# Random seeds
using Random

# Load custom modules
include("$ROOT/scripts/preprocessing/countvars.jl")
include("$ROOT/src/postprocess.jl")
include("$ROOT/src/utilities.jl")
include("$ROOT/scripts/visualization/diagnosticplots.jl")
include("$ROOT/scripts/visualization/paramplots.jl")

# Make custom modules available
using .CountVariables
using .Postprocess
using .CustomUtilityFuns
using .DiagnosticPlots
using .ParamPlots

# Set seed
Random.seed!(42)

# Benchmark model?
benchmark = false

# Load saved chains?
load = false

# Save the result?
save = true
!save && @warn "Samples will NOT be saved automatically."

# If not loading a chain, save results to path below
chainpath = "chains_count_$PRIORSUFFIX.jls"

### MODEL SPECIFICATION ###

lu(x) = length(unique(x))

@model function modelcount(
    r, s, n, PC, y,
    idx_sn, u_n, u_sn, Nv, Ng, Nb;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    α_sxr ~ filldist(Normal(0, σₚ), Ns * Nr)

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, σₚ), Nn, NPC)
    τ_pxn ~ filldist(Exponential(λₚ), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = μ_pxn[u_n, :] .+ τ_pxn[u_n, :] .* z_pxn[u_sn, :]

    σ ~ Exponential(1)

    # Likelihood
    μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
    y ~ MvNormal(μ, σ^2 * I) # Can we LazyArray μ?

    # Generated quantities
    return (; y, α_sxr, β_pxn)
end;

standardise(x) = (x .- mean(x)) ./ std(x)

# Create model
model = modelcount(
    num_region_known,
    num_species_known,
    num_nesting_known,
    PC_known,
    standardise(log.(nbirds)),
    num_species_within_nesting_known,
    unique_nesting_known,
    unique_species_within_nesting_known,
    count_species_by_nesting...,
);

# Benchmark different backends to find out which is fastest
let adbackends = [:forwarddiff, :reversediff, :reversediff_compiled]
    benchmark && benchmark_model(model; check=true, adbackends=adbackends)
end

# Sample from model unless a saved chain should be used
if load
    chain = deserialize("chains/$chainpath")
else
    # Set AD backend to :reversediff and compile with setrdcache(true)
    Turing.setadbackend(:reversediff)
    Turing.setrdcache(true)

    # Configure sampling
    sampler = NUTS(1000, 0.90; max_depth=10)
    nsamples = 2000
    nchains = 4
    ndiscard = 1000

    @info """Sampler: $(string(sampler))
    Samples: $(nsamples)
    Chains: $(nchains)
    Discard: $(ndiscard)
    """

    @info "🚀 Starting sampling: $(Dates.now())"
    chain = sample(model, sampler, MCMCThreads(), nsamples, nchains; discard_initial=ndiscard)

    save && serialize("$ROOT/scripts/models/chains/$chainpath", chain)
    isfile("$ROOT/scripts/models/chains/$chainpath") && @info "💾 Chain saved to '$ROOT/scripts/models/chains/$chainpath'."
end;

θ = generated_quantities(model, chain);

function predictcount(α, β, σ, idx_sn, s, r, X; idx_sr=idx(s, r))
    out = Vector(undef, length(α))
    for i in eachindex(α)
        out[i] = rand.(Normal.(α[i][idx_sr] .+ sum(β[i][idx_sn, :] .* X, dims=2), σ[i]))
    end
    return out
end

function getsamples(θ, sym)
    @chain begin
        map(1:size(θ, 2)) do j
            [θ[i][sym] for i in eachindex(θ[:, j])]
        end
        reduce(vcat, _)
    end
end

α, β = [getsamples(θ, s) for s in [:α_sxr, :β_pxn]];
σ = reduce(hcat, get_params(chain).σ)'

# TODO use [0.75, 0.80, 0.85]
threshold = ppres .> 0.8

countpreds_unknown = @chain begin
    predictcount(
        α,
        β,
        σ,
        num_species_within_nesting_unknown[threshold],
        num_species_unknown[threshold],
        num_region_unknown[threshold],
        PC_unknown[threshold, :],
    )
    reduce(hcat, _)
    Matrix{Float64}(_)
end

countpreds_known = @chain begin
    predictcount(
        α,
        β,
        σ,
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

pred_x = avg_preds_known[num_species_known.==29]
obs_x = log.(nbirds)[num_species_known.==29]
scatter(eachindex(pred_x), pred_x, label="pred")
scatter!(eachindex(pred_x), obs_x, label="obs")

sum(pred_x)
sum(obs_x)

df_countpreds = DataFrame(
    [
        num_atoll_unknown[threshold],
        num_region_unknown[threshold],
        num_species_unknown[threshold],
        avg_preds_unknown
    ],
    [:atoll, :region, :species, :nbirds]
)

CSV.write("$ROOT/data/countpreds_$PRIORSUFFIX.csv", df_countpreds)
@info "Successfully saved predictions to `$ROOT/data/countpreds_$PRIORSUFFIX.csv`."

# Posterior predictive checks
post_preds = let
    model_y_missing = modelcount(
        num_region_known,
        num_species_known,
        num_nesting_known,
        PC_known,
        missing,
        num_species_within_nesting_known,
        unique_nesting_known,
        unique_species_within_nesting_known,
        count_species_by_nesting...,
    )

    generated_quantities(model_y_missing, chain)
end

check_posterior = map(x -> x.y, vec(post_preds))

histogram(log.(nbirds), normalize=true, alpha=0.5)
histogram!(vcat(check_posterior...), normalize=true, c=:white, lw=3)
histogram!(vcat(check_posterior...), normalize=true, c=2, lw=1.5, legend=:none)

# Prior predictive checks ... not yet
model_missing = modelcount(
    [missing for _ in eachindex(num_region_known)],
    [missing for _ in eachindex(num_species_known)],
    [missing for _ in eachindex(num_nesting_known)],
    reshape([missing for _ in eachindex(PC_known)], size(PC_known)),
    [missing for _ in eachindex(nbirds)],
    num_species_within_nesting_known,
    unique_nesting_known,
    unique_species_within_nesting_known,
    count_species_by_nesting...,
    Nr=lu(num_region_known),
    Ns=lu(num_species_known),
    Nn=lu(num_nesting_known),
    NPC=size(PC_known, 2),
    idx_sr=idx(num_species_known, num_region_known)
)

# Prior predictive checks
prior_preds = let
    prior_chain = sample(model, Prior(), 5000)
    predict(model_missing, prior_chain)
end;

density(get_params(prior_chain).τ_sxr)

check_prior = map(x -> x.y, vec(prior_preds))

histogram(log.(nbirds), normalize=true, lw=1.5, bins=20)
histogram!(vcat(check_prior...), normalize=true, c=2, lw=1.5, legend=:none, bins=20)
