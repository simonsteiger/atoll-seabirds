# Pathing
const ROOT = dirname(Base.active_project())

# Probabilistic programming
using Turing, TuringBenchmarking, ReverseDiff, ParetoSmooth
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
include("$ROOT/src/presence.jl")
include("$ROOT/src/utilities.jl")

# Make custom modules available
using .PresenceVariables
using .CustomUtilityFuns

# Set seed
Random.seed!(42)

# Benchmark model?
benchmark = false

# Run the sampler?
run = isempty(ARGS) ? false : ARGS[1] == "true"

### MODEL SPECIFICATION ###

# Prior settings can be set from command line
σₚ, θₚ = 1, 0.1

if !run
    @info "Loading chain, no model fit."
elseif isempty(ARGS) || ARGS[2] == "default"
    @info "Fitting model with default priors: σₚ=$σₚ, θₚ=$θₚ."
elseif all(ARGS[2] .!= ["narrow", "wide"])
    throw("Unknown prior setting: '$(ARGS[2])'. Pass nothing or one of 'default', 'narrow', 'wide'.")
else
    σₚ, θₚ = ARGS[2] == "wide" ? [σₚ, θₚ] .* 3 : [σₚ, θₚ] .* 1 / 3
    @info "Fitting model with $(ARGS[2]) priors: σₚ=$(round(σₚ, digits=2)), θₚ=$(round(θₚ, digits=2))."
end

PRIORSUFFIX = isempty(ARGS) ? "default" : ARGS[2]

# Path to read or write chain from / to
chainpath = "presence_$PRIORSUFFIX.jls"

lu(x) = length(unique(x))

function peaceful_generated_quantities(m, c)
    chains_params = Turing.MCMCChains.get_sections(c, :parameters)
    return generated_quantities(m, chains_params)
end

@model function modelpresence(
    r, s, n, PC, y,
    idx_sn, u_n, u_sn, Nv, Ng, Nb;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    α_sxr ~ filldist(Normal(0, σₚ), Ns * Nr)

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, σₚ), Nn, NPC)
    τ_pxn ~ filldist(InverseGamma(3, θₚ), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = @. μ_pxn[u_n, :] + τ_pxn[u_n, :] * z_pxn[u_sn, :]

    # Likelihood
    v = α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2)
    y .~ BernoulliLogit.(v)

    # Generated quantities
    return (; y, α_sxr, β_pxn)
end;

# Create model
model = modelpresence(
    num_region,
    num_species,
    num_nesting,
    PC,
    presence,
    num_species_within_nesting,
    unique_nesting,
    unique_species_within_nesting,
    count_species_by_nesting...,
);

# Prior predictive checks
check_prior = sample(model, Prior(), 5000)

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
    nsamples = 5000
    nthreads = 4
    ndiscard = 500

    @info """Sampler: $(string(sampler))
    Samples: $(nsamples)
    Threads: $(nthreads)
    Discard: $(ndiscard)
    """

    @info "🚀 Starting sampling: $(Dates.now())"
    chain = sample(model, sampler, MCMCThreads(), nsamples, nthreads; discard_initial=ndiscard)

    serialize("$ROOT/results/chains/$chainpath", chain)
    @info "💾 Chain saved to `$ROOT/results/chains/$chainpath`."
end;

θ = peaceful_generated_quantities(model, chain)

function predictpresence(α, β, idx_sn, s, r, X; idx_sr=idx(s, r))
    [rand.(BernoulliLogit.(α[i][idx_sr] .+ sum(β[i][idx_sn, :] .* X, dims=2))) for i in eachindex(α)]
end

α = [θ[i].α_sxr for i in eachindex(θ[:, 1])]
β = [θ[i].β_pxn for i in eachindex(θ[:, 1])]

nnu_long = [fill(num_species_within_nesting_unknown, length(num_region_unknown))...;]
nru_long = [fill.(num_region_unknown, length(num_nesting_unknown))...;]
nsu_long = [fill(num_species_unknown, length(num_region_unknown))...;]
PCu_long = reduce(vcat, [permutedims(hcat(fill(s, length(num_nesting_unknown))...)) for s in eachslice(PC_unknown, dims=1)])

preds = Matrix{Float64}(reduce(hcat, vec(predictpresence(α, β, nnu_long, nsu_long, nru_long, PCu_long));))

pct_preds = vec(mean(preds, dims=2))

ssu_long = [fill(PresenceVariables.str_species_unknown, length(num_region_unknown))...;]

sau_long = [fill.(PresenceVariables.str_atoll_unknown, length(num_nesting_unknown))...;]

df_preds = @chain begin
    DataFrame([sau_long, ssu_long, pct_preds], [:atoll, :species, :percent])
    unstack(_, :species, :percent)
end

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