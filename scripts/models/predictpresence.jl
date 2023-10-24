# Full model to predict species presence

### SETUP ###

cd("scripts/models")

# Probabilistic programming
using Turing, TuringBenchmarking, LazyArrays
# Statistics
using StatsFuns, LinearAlgebra
# Working with tabular data
using Chain, DataFrames
# Plotting
using StatsPlots
# Saving results
using Serialization, CSV

# Load custom modules
include("../preprocessing/preprocess.jl")
include("../../src/postprocess.jl")
include("../../src/utilities.jl")
include("../visualization/diagnosticplots.jl")
include("../visualization/paramplots.jl")

# Add custom modules
using .Preprocess
using .Postprocess
using .CustomUtilityFuns
using .DiagnosticPlots
using .ParamPlots

# Benchmark model?
benchmark = false

# Load a saved chain?
load = false

# If not loading a chain, save results to path below
savetofile = "predictpresence_new"

# If loading a chain, which file should be loaded?
loadfrompath = "chains/predictpresence_new.jls"

### MODEL SPECIFICATION ###

lu(x) = length(unique(x))
idx(i, j) = i .+ (j .- 1) * maximum(i)

@model function modelpresence(r, s, n, PC, y;       # Main inputs
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2),  # Number of groups
    idx_sr=idx(s, r), idx_PC=[1:NPC;])              # Vectors for indexing

    # Priors for species × region
    μ_sxr ~ Normal()
    τ_sxr ~ Exponential(2)
    z_sxr ~ filldist(Normal(), Ns * Nr)
    α_sxr = μ_sxr .+ τ_sxr .* getindex.((z_sxr,), idx_sr)

    # Priors for burrow nesters × PCs
    μ_pxb ~ Normal()
    τ_pxb ~ Exponential(2)
    z_pxb ~ filldist(Normal(), NPC)
    β_pxb = μ_pxb .+ τ_pxb .* getindex.((z_pxb,), idx_PC)

    # Priors for ground nesters × PCs
    μ_pxg ~ Normal()
    τ_pxg ~ Exponential(2)
    z_pxg ~ filldist(Normal(), NPC)
    β_pxg = μ_pxg .+ τ_pxg .* getindex.((z_pxg,), idx_PC)

    # Priors for tree nesters × PCs
    μ_pxv ~ Normal()
    τ_pxv ~ Exponential(2)
    z_pxv ~ filldist(Normal(), NPC)
    β_pxv = μ_pxv .+ τ_pxv .* getindex.((z_pxv,), idx_PC)

    # Convert to matrix for vectorization
    β_pxn = reshape([β_pxb; β_pxg; β_pxv], (NPC, Nn))

    # Likelihood
    v = α_sxr[idx_sr] .+ β_pxn[:, n]' .* PC
    y .~ BernoulliLogit.(v)

    return (α_sxr=α_sxr, β_pxn=β_pxn)
end;

model = modelpresence(
    num_region,
    num_species,
    num_nesting,
    PC,
    presence
);

### SAMPLING ###

# Benchmark different backends to find out which is fastest
if benchmark
    adbackends = [:forwarddiff, :reversediff, :reversediff_compiled]
    benchmark_model(
        model;
        check=true, # Assert that model is correctly specified
        adbackends=adbackends
    )
end

# ------- BENCHMARK RESULTS ------- #
# 18.817 ms ReverseDiff             #
# 28.948 ms ForwardDiff             #
#  6.608 ms ReverseDiff[compiled]   #
# ??.??? ms Zygote "never" finished #
# --------------------------------- #

# Set AD backend to :reversediff and compile with setrdcache(true)
Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

# Sample from model
if load
    chain = deserialize(loadfrompath)
else
    sampler = NUTS(1000, 0.65; max_depth=6)
    chain = sample(
        model,
        sampler,
        MCMCThreads(),
        2_000,
        3;
        discard_initial=1000
    )
    serialize("chains/$savetofile.jls", chain)
end

# Check acceptance rate
plot_acceptance_rate(chain)
# Check rhat
plot_rhat(chain)

θs = ["θ$i$j" for i in 1:6, j in 1:3]

params = extractparams(chain, ["λ", θs...])

preds = Dict()

for i in eachindex(num_species_unknown)
    preds[Preprocess.str_species_unknown[i]] =
        prediction(
            "presence",
            params,
            num_species_unknown[i],
            num_species_within_nesting_unknown[i],
            num_nesting_unknown[i],
            num_region_unknown,
            num_atoll_unknown,
            PC_unknown
        )
end

for k in keys(preds)
    preds[k] = mean.(preds[k])
end

df = DataFrame([preds...])
insertcols!(df, 1, :atoll => Preprocess.str_atoll_unknown)

CSV.write("../../data/$savetofile.csv", df)

function specieshisto(df, idx)
    p = histogram(df[:, idx], bins=10)
    xlims!(0, 1)
    title!(names(df)[idx])
    return p
end

specieshisto(df, 37)

# TODO
# Add human population? PC or directly?
# Try TDist
# Prior predictive simulation - how sensitive is our model to variation in the predictors? Shouldn't be suuuper sensitive
# LOO - how sensitive is our model to single data points in the known data
