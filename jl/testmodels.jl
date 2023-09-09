using Turing, StatsFuns, StatsPlots

include("preprocess.jl")
include("tune.jl")

using .Preprocess

# random intercept only
@model function random_intercept(atoll, presence)
    a = length(unique(atoll))
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), a)

    for i in eachindex(presence)
        v = logistic(α[atoll[i]])
        presence[i] ~ Bernoulli(v)
    end
end;

m1 = random_intercept(all_atoll, all_presence)
chain1 = sample(m1, HMC(0.05, 10), 15_000, discard_initial=5000)

# species and atoll intercept
@model function species_intercept(atoll, species, presence)
    a = length(unique(atoll))
    s = length(unique(species))
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), a)
    β ~ filldist(Normal(0, 1), s)

    for i in eachindex(presence)
        v = logistic(α[atoll[i]] + β[species[i]])
        presence[i] ~ Bernoulli(v)
    end
end;

m2 = species_intercept(all_atoll, all_species, all_presence)
chain2 = sample(m2, HMC(0.05, 10), 8000, discard_initial=2000)

@model function pc_species_intercept(atoll, species, pc, presence)
    # Number of groups per predictor
    a = length(unique(atoll))
    s = length(unique(species))

    # Priors for atolls
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), a)
    # Priors for species effect
    β ~ filldist(Normal(0, 1), s)
    # Priors for PC
    θ1 ~ filldist(Normal(0, 1), s)
    θ2 ~ filldist(Normal(0, 1), s)
    θ3 ~ filldist(Normal(0, 1), s)
    θ4 ~ filldist(Normal(0, 1), s)
    θ5 ~ filldist(Normal(0, 1), s)
    θ6 ~ filldist(Normal(0, 1), s)

    for i in eachindex(presence)
        v = logistic(
            α[atoll[i]] +
            β[species[i]] +
            θ1[species[i]] * pc[i, 1] +
            θ2[species[i]] * pc[i, 2] +
            θ3[species[i]] * pc[i, 3] +
            θ4[species[i]] * pc[i, 4] +
            θ5[species[i]] * pc[i, 5] +
            θ6[species[i]] * pc[i, 6]
        )
        presence[i] ~ Bernoulli(v)
    end
end;

m3 = pc_species_intercept(all_atoll, all_species, all_pc, all_presence)
chain3 = sample(m3, HMC(0.05, 10), 8000, discard_initial=2000)
