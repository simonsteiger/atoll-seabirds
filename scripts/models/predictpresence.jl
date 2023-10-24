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
using Serialization

# Load custom modules
include("../preprocessing/preprocess.jl")
include("../../src/postprocess.jl")
include("../../src/utilities.jl")
include("../visualization/diagnosticplots.jl")
include("../visualization/paramplots.jl")

# Add names to environment
using .Preprocess
using .Postprocess
using .CustomUtilityFuns
using .DiagnosticPlots
using .ParamPlots

# Benchmark model?
benchmark = false

# Load a saved chain?
load = true

# Which chain should be loaded?
filepath = "chains/predictpresence.jls"

### MODEL SPECIFICATION ###

@model function modelpresence(PC, r, s, n, s_in_n, Nv, Ng, Nb, y)
    # Number of groups per predictor
    Nr = length(unique(r))
    Ns = length(unique(s))
    Nn = length(unique(n))

    # Priors for region
    λ̄ ~ Normal()
    κ ~ Exponential(1)
    λ ~ filldist(Normal(λ̄, κ), Nr, Ns)

    # PC1 per nesting type
    θ̄1 ~ filldist(Normal(), Nn) # Distribution of distributions per nesting type
    τ1 ~ filldist(Exponential(1), Nn)
    θ11 ~ filldist(Normal(θ̄1[1], τ1[1]), Nb) # Vector of distributions of ground nesters
    θ12 ~ filldist(Normal(θ̄1[2], τ1[2]), Ng) # Vector of distributions of ground nesters
    θ13 ~ filldist(Normal(θ̄1[3], τ1[3]), Nv) # Vector of distributions of ground nesters
    θ1 = [θ11, θ12, θ13]

    # PC2 per nesting type
    θ̄2 ~ filldist(Normal(), Nn) # Distribution of distributions per nesting type
    τ2 ~ filldist(Exponential(1), Nn)
    θ21 ~ filldist(Normal(θ̄2[1], τ2[1]), Nb) # Vector of distributions of ground nesters
    θ22 ~ filldist(Normal(θ̄2[2], τ2[2]), Ng) # Vector of distributions of tree nesters
    θ23 ~ filldist(Normal(θ̄2[3], τ2[3]), Nv) # Vector of distributions of burrow nesters
    θ2 = [θ21, θ22, θ23]

    # PC3 per nesting type
    θ̄3 ~ filldist(Normal(), Nn) # Distribution of distributions per nesting type
    τ3 ~ filldist(Exponential(1), Nn)
    θ31 ~ filldist(Normal(θ̄3[1], τ3[1]), Nb) # Vector of distributions of ground nesters
    θ32 ~ filldist(Normal(θ̄3[2], τ3[2]), Ng) # Vector of distributions of tree nesters
    θ33 ~ filldist(Normal(θ̄3[3], τ3[3]), Nv) # Vector of distributions of burrow nesters
    θ3 = [θ31, θ32, θ33]

    # PC4 per nesting type
    θ̄4 ~ filldist(Normal(), Nn) # Distribution of distributions per nesting type
    τ4 ~ filldist(Exponential(1), Nn)
    θ41 ~ filldist(Normal(θ̄4[1], τ4[1]), Nb) # Vector of distributions of ground nesters
    θ42 ~ filldist(Normal(θ̄4[2], τ4[2]), Ng) # Vector of distributions of tree nesters
    θ43 ~ filldist(Normal(θ̄4[3], τ4[3]), Nv) # Vector of distributions of burrow nesters
    θ4 = [θ41, θ42, θ43]

    # PC5 per nesting type
    θ̄5 ~ filldist(Normal(), Nn) # Distribution of distributions per nesting type
    τ5 ~ filldist(Exponential(1), Nn)
    θ51 ~ filldist(Normal(θ̄5[1], τ5[1]), Nb) # Vector of distributions of ground nesters
    θ52 ~ filldist(Normal(θ̄5[2], τ5[2]), Ng) # Vector of distributions of tree nesters
    θ53 ~ filldist(Normal(θ̄5[3], τ5[3]), Nv) # Vector of distributions of burrow nesters
    θ5 = [θ51, θ52, θ53]

    # PC6 per nesting type
    θ̄6 ~ filldist(Normal(), Nn) # Distribution of distributions per nesting type
    τ6 ~ filldist(Exponential(1), Nn)
    θ61 ~ filldist(Normal(θ̄6[1], τ6[1]), Nb) # Vector of distributions of ground nesters
    θ62 ~ filldist(Normal(θ̄6[2], τ6[2]), Ng) # Vector of distributions of tree nesters
    θ63 ~ filldist(Normal(θ̄6[3], τ6[3]), Nv) # Vector of distributions of burrow nesters
    θ6 = [θ61, θ62, θ63]

    for i in eachindex(y)
        p = logistic(
            λ[r[i], s[i]] +
            θ1[n[i]][s_in_n[i]] * PC[i, 1] +
            θ2[n[i]][s_in_n[i]] * PC[i, 2] +
            θ3[n[i]][s_in_n[i]] * PC[i, 3] +
            θ4[n[i]][s_in_n[i]] * PC[i, 4] +
            θ5[n[i]][s_in_n[i]] * PC[i, 5] +
            θ6[n[i]][s_in_n[i]] * PC[i, 6]
        )
        y[i] ~ Bernoulli(p)
    end
end;

model = modelpresence(
    PC,
    num_region,
    num_species,
    num_nesting,
    num_species_within_nesting,
    count_species_by_nesting...,
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
    chain = deserialize(filepath);
else
    chain = sample(
        model,
        HMC(0.025, 10), # tuned to acceptance_rate ≈ 0.65, see https://pythonhosted.org/pyhmc/tuning.html
        MCMCThreads(),
        40_000,         # number of samples
        4;              # number of chains
        discard_initial=2000
    )
    serialize("chains/predictpresence.jls", chain)
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

using CSV

# CSV.write("../../data/predictpresence.csv", df)

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
