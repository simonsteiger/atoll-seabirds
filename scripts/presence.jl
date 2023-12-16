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
using Serialization, CSV, Dates, Markdown
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

# Benchmark model?
benchmark = false

# Run the sampler?
sample_model = isempty(ARGS) ? false : ARGS[1] == "run"

if ARGS[2] == "S" # simplest model
    model = m1p
    inputs = (; num_region, num_species)
elseif isempty(ARGS) || ARGS[2] == "M" # most complex model
    model = mMp
    # What does the auto named tuple notation (; ...) do with splatted vectors?
    inputs = (num_region, num_species, num_nesting, PC, num_species_within_nesting, unique_nesting, unique_species_within_nesting, count_species_by_nesting...,)
elseif isempty(ARGS) || ARGS[2] == "L" # large model
    model = mLp
    # What does the auto named tuple notation (; ...) do with splatted vectors?
    inputs = (num_region, num_species, num_nesting, PC, num_species_within_nesting, unique_nesting, unique_species_within_nesting, count_species_by_nesting...,)
elseif ARGS[2] == "XL"
    model = mXLp
    # What does the auto named tuple notation (; ...) do with splatted vectors?
    inputs = (num_atoll, num_region, num_species, num_nesting, PC, num_species_within_nesting, unique_nesting, unique_species_within_nesting, count_species_by_nesting...,)
end

# Prior settings can be set from command line
# σₚ, θₚ = 1, 1

# if !sample_model
#     @info "Loading chain, no model fit."
# elseif isempty(ARGS) || ARGS[2] == "default"
#     @info "Fitting model with default priors: σₚ=$σₚ, θₚ=$θₚ."
# elseif all(ARGS[2] .!= ["narrow", "wide"])
#     throw("Unknown prior setting: '$(ARGS[2])'. Pass nothing or one of 'default', 'narrow', 'wide'.")
# else
#     σₚ, θₚ = ARGS[2] == "wide" ? [σₚ, θₚ] .* 3 : [σₚ, θₚ] .* 1 / 3
#     @info "Fitting model with $(ARGS[2]) priors: σₚ=$(round(σₚ, digits=2)), θₚ=$(round(θₚ, digits=2))."
# end

# PRIORSUFFIX = isempty(ARGS) ? "default" : ARGS[2]

# Path to read or write chain from / to
chainpath = "presence.jls"

function peaceful_generated_quantities(m, c)
    chains_params = Turing.MCMCChains.get_sections(c, :parameters)
    return generated_quantities(m, chains_params)
end

# --- PRIOR PREDICTIVE CHECK --- #

chain_prior = sample(model(inputs..., presence), Prior(), 1000);

θ_prior = peaceful_generated_quantities(model(inputs..., presence), chain_prior);

let k = keys(θ_prior[1])
    params = vec(getsamples(θ_prior, k...))
    priorsamples = reduce(vcat, simulate(params, inputs...))
    density(priorsamples, fillrange=0, fillalpha=0.2, legend=false)
end

# Benchmark different backends to find out which is fastest
if benchmark
    backends = [
           Turing.Essential.ForwardDiffAD{0}(),
           Turing.Essential.ReverseDiffAD{false}(),
           Turing.Essential.ReverseDiffAD{true}()
       ];
    
    TuringBenchmarking.run(TuringBenchmarking.make_turing_suite(model(inputs..., presence), adbackends=backends);)
end

# Sample from model unless a saved chain should be used
if !sample_model
    # chain = deserialize("$ROOT/results/chains/$chainpath")
else
    # Set AD backend to :reversediff and compile with setrdcache(true)
    Turing.setadbackend(:reversediff)
    Turing.setrdcache(true)

    # Configure sampling
    sampler = NUTS(1000, 0.95; max_depth=10)
    nsamples = 2000
    nthreads = 4
    ndiscard = 200

    @info """Sampler: $(string(sampler))
    Samples: $(nsamples)
    Threads: $(nthreads)
    Discard: $(ndiscard)
    """

    @info "🚀 Starting sampling: $(Dates.now())"
    chain = sample(model(inputs..., presence), sampler, MCMCThreads(), nsamples, nthreads; discard_initial=ndiscard)

    # serialize("$ROOT/results/chains/$chainpath", chain)
    # @info "💾 Chain saved to `$ROOT/results/chains/$chainpath`."
end;

# --- POSTERIOR PREDICTIVE CHECK --- #

θ_post = peaceful_generated_quantities(model(inputs..., presence), chain)

postpcplots = let k = keys(θ_post[1])
    params = vec(getsamples(θ_post, k...))
    postsamples = reduce(hcat, simulate(params, inputs...))

    map(enumerate(unique(num_species))) do (index, species)
        pred_x = vec(postsamples[num_species.==species, :])
        obs_x = presence[num_species.==species]
        histogram(obs_x, normalize=true, alpha=0.5, lc=:transparent, bins=10, label="O")
        density!(pred_x, fillrange=0, fillalpha=0.2, normalize=true, alpha=0.5, bins=10, lc=:transparent, yticks=:none, label="P")
        xticks!(0:0.5:1, string.(0:0.5:1)); 
        title!(unique(str_species)[index], titlefontsize=8)
        vline!([0.8], c=:black, ls=:dash, label=:none)
    end
end

plot(postpcplots..., titlefontsize=9, size=(800,1200), layout=(8,5))
savefig("results/svg/postpreds_L.svg")


params_post = []

for sym in [:α_sxr, :β_n, :β_pxn]
    res = haskey(θ_post[1], sym) && vec(map(x -> getfield(x, sym), θ_post))
    res != false && push!(params_post, res)
end

postpc = @chain begin
    simulate(params_post..., inputs...)
    reduce(hcat, _)
    Matrix{Float64}(_)
    mean(_, dims=2)
end

hist_species_postpc = map(enumerate(unique(num_species))) do (index, species)
    pred_x = postpc[num_species.==species]
    obs_x = presence[num_species.==species]
    histogram(obs_x, normalize=true, alpha=0.5, lc=:transparent, bins=10, label="O")
    histogram!(pred_x, normalize=true, alpha=0.5, bins=10, lc=:transparent, yticks=:none, label="P")
    xticks!(0:0.5:1, string.(0:0.5:1)); 
    title!(unique(str_species)[index], titlefontsize=8)
end

plot(hist_species_postpc..., size=(800,1200), layout=(8,5))

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