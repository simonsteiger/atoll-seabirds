# Full model to predict species presence

### SETUP ###

cd("scripts/models")

# Probabilistic programming
using Turing, TuringBenchmarking
# Statistics
using StatsFuns, LinearAlgebra
# Working with tabular data
using Chain, DataFrames
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
benchmark = true

# Load a saved chain?
load = false

# If not loading a chain, save results to path below
savetofile = "predictpresence_T"

# If loading a chain, which file should be loaded?
loadfrompath = "chains/predictpresence.jls"

### MODEL SPECIFICATION ###

lu(x) = length(unique(x))

@model function modelpresence(r, s, y; Nr=lu(r), Ns=lu(s), k=maximum(s))
    # Priors for region
    θ_reg ~ filldist(TDist(3), Nr)

    # Priors for species
    θ_spe ~ filldist(TDist(3), Ns)

    # Priors for region x species
    τ_rxs ~ Exponential(0.5)
    θ_rxs ~ filldist(Normal(0, τ_rxs), Nr * Ns)

    for i in eachindex(y)
        p = logistic(θ_reg[r[i]] + θ_spe[s[i]] + θ_rxs[s[i]+(r[i]-1)*k])
        y[i] ~ Bernoulli(p)
    end
end;

model = modelpresence(
    num_region,
    num_species,
    presence
);

### SAMPLING ###

# Benchmark different backends to find out which is fastest
if benchmark
    adbackends = [:forwarddiff, :reversediff, :reversediff_compiled]
    benchmark_model(model; check=true, adbackends=adbackends)
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
    chain = deserialize(loadfrompath);
else
    chain = sample(
        model,
        HMC(0.025, 10), # tuned to acceptance_rate ≈ 0.65, see https://pythonhosted.org/pyhmc/tuning.html
        #MCMCThreads(),
        5_000,         # number of samples
        #1;              # number of chains
        discard_initial=2000
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
