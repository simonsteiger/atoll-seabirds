# Full model to predict species presence

using Turing, StatsFuns, StatsPlots, Chain, DataFrames

include("../preprocessing/preprocess.jl")
include("../../src/postprocess.jl")

using .Preprocess
using .Postprocess
using ReverseDiff
using StatsBase
using LinearAlgebra

# For saving stuff
using Serialization, Dates

const TODAY = Dates.format(now(), "yyyy-mm-dd")

Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

@model function full_model(PC, r, s, n, s_in_n, Nv, Ng, Nb, y)
    # Number of groups per predictor
    Nr = length(unique(r))
    Ns = length(unique(s))
    Nn = length(unique(n))

    # Priors for region
    λ̄ ~ Normal()
    κ ~ Exponential(1)
    λ ~ filldist(Normal(λ̄[1], κ[1]), Nr, Ns)

    # PC1 per nesting type
    θ̄1 ~ filldist(Normal(), Nn) # Distribution of distributions of ground nesters
    τ1 ~ filldist(Exponential(1), Nn)
    θ11 ~ filldist(Normal(θ̄1[1], τ1[1]), Nb) # Vector of distributions of ground nesters
    θ12 ~ filldist(Normal(θ̄1[2], τ1[2]), Ng) # Vector of distributions of tree nesters
    θ13 ~ filldist(Normal(θ̄1[3], τ1[3]), Nv) # Vector of distributions of burrow nesters
    θ1 = [θ11, θ12, θ13]

    # PC2 per nesting type
    θ̄2 ~ filldist(Normal(), Nn) # Distribution of distributions of ground nesters
    τ2 ~ filldist(Exponential(1), Nn)
    θ21 ~ filldist(Normal(θ̄2[1], τ2[1]), Nb) # Vector of distributions of ground nesters
    θ22 ~ filldist(Normal(θ̄2[2], τ2[2]), Ng) # Vector of distributions of tree nesters
    θ23 ~ filldist(Normal(θ̄2[3], τ2[3]), Nv) # Vector of distributions of burrow nesters
    θ2 = [θ21, θ22, θ23]

    # PC3 per nesting type
    θ̄3 ~ filldist(Normal(), Nn) # Distribution of distributions of ground nesters
    τ3 ~ filldist(Exponential(1), Nn)
    θ31 ~ filldist(Normal(θ̄3[1], τ3[1]), Nb) # Vector of distributions of ground nesters
    θ32 ~ filldist(Normal(θ̄3[2], τ3[2]), Ng) # Vector of distributions of tree nesters
    θ33 ~ filldist(Normal(θ̄3[3], τ3[3]), Nv) # Vector of distributions of burrow nesters
    θ3 = [θ31, θ32, θ33]

    # PC4 per nesting type
    θ̄4 ~ filldist(Normal(), Nn) # Distribution of distributions of ground nesters
    τ4 ~ filldist(Exponential(1), Nn)
    θ41 ~ filldist(Normal(θ̄4[1], τ4[1]), Nb) # Vector of distributions of ground nesters
    θ42 ~ filldist(Normal(θ̄4[2], τ4[2]), Ng) # Vector of distributions of tree nesters
    θ43 ~ filldist(Normal(θ̄4[3], τ4[3]), Nv) # Vector of distributions of burrow nesters
    θ4 = [θ41, θ42, θ43]

    # PC5 per nesting type
    θ̄5 ~ filldist(Normal(), Nn) # Distribution of distributions of ground nesters
    τ5 ~ filldist(Exponential(1), Nn)
    θ51 ~ filldist(Normal(θ̄5[1], τ5[1]), Nb) # Vector of distributions of ground nesters
    θ52 ~ filldist(Normal(θ̄5[2], τ5[2]), Ng) # Vector of distributions of tree nesters
    θ53 ~ filldist(Normal(θ̄5[3], τ5[3]), Nv) # Vector of distributions of burrow nesters
    θ5 = [θ51, θ52, θ53]

    # PC6 per nesting type
    θ̄6 ~ filldist(Normal(), Nn) # Distribution of distributions of ground nesters
    τ6 ~ filldist(Exponential(1), Nn)
    θ61 ~ filldist(Normal(θ̄6[1], τ6[1]), Nb) # Vector of distributions of ground nesters
    θ62 ~ filldist(Normal(θ̄6[2], τ6[2]), Ng) # Vector of distributions of tree nesters
    θ63 ~ filldist(Normal(θ̄6[3], τ6[3]), Nv) # Vector of distributions of burrow nesters
    θ6 = [θ61, θ62, θ63]

    for i in eachindex(presence)
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

model = full_model(
    PC,
    num_region,
    num_species,
    num_nesting,
    num_species_within_nesting,
    count_species_by_nesting...,
    presence
);

chain3 = sample(model, HMC(0.01, 10), MCMCThreads(), 30_000, 4) # NUTS ETA was 19:55h

serialize("$(TODAY)_chain_joinregion.jls", chain)

params = Postprocess.extractparams(chain, Postprocess.targets)

λ = @chain chain3 begin
    group(_, "λ")
    mean(_)
    _[:, 2]
    reshape(_, 4, 37)
end

thetas = [mean(group(chain3, "θ$i$j"))[:, 2] for i in 1:6, j in 1:3]
θ11, θ12, θ13 = [thetas[1, i] for i in 1:3]
θ21, θ22, θ23 = [thetas[1, i] for i in 1:3]
θ31, θ32, θ33 = [thetas[1, i] for i in 1:3]
θ41, θ42, θ43 = [thetas[1, i] for i in 1:3]
θ51, θ52, θ53 = [thetas[1, i] for i in 1:3]
θ61, θ62, θ63 = [thetas[1, i] for i in 1:3]

params = (
    λ,
    θ1=[θ11, θ12, θ13],
    θ2=[θ21, θ22, θ23],
    θ3=[θ31, θ32, θ33],
    θ4=[θ41, θ42, θ43],
    θ5=[θ51, θ52, θ53],
    θ6=[θ61, θ62, θ63],
)

preds = Dict{Any,Vector{Bool}}()

for i in eachindex(Preprocess.num_species_unknown)
    preds[Preprocess.str_species_unknown[i]] =
        Postprocess.prediction(
            "presence",
            params,
            Preprocess.num_species_unknown[i],
            Preprocess.num_species_within_nesting_unknown[i],
            Preprocess.num_nesting_unknown[i],
            Preprocess.num_region_unknown,
            Preprocess.num_atoll_unknown,
            Preprocess.PC_unknown
        )
end

[println("$k: $(round(sum(preds[k])/84*100, digits=1))") for k in keys(preds)]

# TODO
# Use all samples for prediction
# Plot all atolls color coded known true/false and prediction certainty
