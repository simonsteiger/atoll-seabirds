using Turing, StatsFuns, StatsPlots, Chain, DataFrames

include("preprocess.jl")
include("tune.jl")

using .Preprocess
using ReverseDiff

Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

function printchain(chain)
    x = describe(chain)[1]
    show(DataFrame(x), allrows=true)
end

# # random intercept only
# @model function random_intercept(atoll, presence)
#     a = length(unique(atoll))
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
# 
#     for i in eachindex(presence)
#         v = logistic(α[atoll[i]])
#         presence[i] ~ Bernoulli(v)
#     end
# end;
# 
# m1 = random_intercept(all_atoll, all_presence)
# #chain1 = sample(m1, HMC(0.05, 10), 15_000, discard_initial=5000)
# 
# # species and atoll intercept
# @model function species_intercept(atoll, species, presence)
#     a = length(unique(atoll))
#     s = length(unique(species))
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
#     β ~ filldist(Normal(0, 1), s)
# 
#     for i in eachindex(presence)
#         v = logistic(α[atoll[i]] + β[species[i]])
#         presence[i] ~ Bernoulli(v)
#     end
# end;
# 
# m2 = species_intercept(all_atoll, all_species, all_presence)
# #chain2 = sample(m2, HMC(0.05, 10), 8000, discard_initial=2000)
# 
# @model function pc_species_intercept(atoll, species, pc, presence)
#     # Number of groups per predictor
#     a = length(unique(atoll))
#     s = length(unique(species))
# 
#     # Priors for atolls
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
#     # Priors for species effect
#     β ~ filldist(Normal(0, 1), s)
#     # Priors for PC
#     θ1 ~ filldist(Normal(0, 1), s)
#     θ2 ~ filldist(Normal(0, 1), s)
#     θ3 ~ filldist(Normal(0, 1), s)
#     θ4 ~ filldist(Normal(0, 1), s)
#     θ5 ~ filldist(Normal(0, 1), s)
#     θ6 ~ filldist(Normal(0, 1), s)
# 
#     for i in eachindex(presence)
#         v = logistic(
#             α[atoll[i]] +
#             β[species[i]] +
#             θ1[species[i]] * pc[i, 1] +
#             θ2[species[i]] * pc[i, 2] +
#             θ3[species[i]] * pc[i, 3] +
#             θ4[species[i]] * pc[i, 4] +
#             θ5[species[i]] * pc[i, 5] +
#             θ6[species[i]] * pc[i, 6]
#         )
#         presence[i] ~ Bernoulli(v)
#     end
# end;

using StatsBase

num_species = @chain Preprocess.envs_known begin
    groupby(_, :nestingtype)
    DataFrames.transform(_, :species => denserank => :num_species)
    DataFrames.transform(_, :num_species => ByRow(x -> Int64(x)) => identity)
    getproperty(_, :num_species)
end

vgb = @chain Preprocess.envs_known begin
    unique(_, :species)
    groupby(_, :nestingtype)
    combine(_, nrow)
    getproperty(_, :nrow)
end

# get species effect back in
# dont throw out the albatrosses (dont exclude)
# all 6 pcs

@model function pc_species_nesting(atoll, species, nesting, v, g, b, pc, presence)
    # Number of groups per predictor
    a = length(unique(atoll))
    #s = length(unique(species))

    # Priors for atolls
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), a)

    # Priors for species effect
    # β ~ filldist(Normal(0, 1), s)

    ω̄11 ~ Normal() # Distribution of distributions of ground nesters
    ω̄12 ~ Normal()
    ω̄13 ~ Normal()
    κ11 ~ Exponential(1)
    κ12 ~ Exponential(1)
    κ13 ~ Exponential(1)
    ω11 ~ filldist(Normal(ω̄11, κ11), b) # Vector of distributions of ground nesters
    ω12 ~ filldist(Normal(ω̄12, κ12), g) # Vector of distributions of tree nesters
    ω13 ~ filldist(Normal(ω̄13, κ13), v) # Vector of distributions of burrow nesters
    ω1 = [ω11, ω12, ω13]

    ω̄21 ~ Normal() # Distribution of distributions of ground nesters
    ω̄22 ~ Normal()
    ω̄23 ~ Normal()
    κ21 ~ Exponential(1)
    κ22 ~ Exponential(1)
    κ23 ~ Exponential(1)
    ω21 ~ filldist(Normal(ω̄21, κ21), b) # Vector of distributions of ground nesters
    ω22 ~ filldist(Normal(ω̄22, κ22), g) # Vector of distributions of tree nesters
    ω23 ~ filldist(Normal(ω̄23, κ23), v) # Vector of distributions of burrow nesters
    ω2 = [ω21, ω22, ω23]

    for i in eachindex(presence)
        v = logistic(
            α[atoll[i]] +
            #β[species[i]] +
            ω1[nesting[i]][species[i]] * pc[i, 1] +
            ω2[nesting[i]][species[i]] * pc[i, 2] #+
            #θ3[nesting[i], species[i]] * pc[i, 3] +
            #θ4[nesting[i], species[i]] * pc[i, 4] +
            #θ5[nesting[i], species[i]] * pc[i, 5] +
            #θ6[nesting[i], species[i]] * pc[i, 6]
        )
        presence[i] ~ Bernoulli(v)
    end
end;

m3 = pc_species_nesting(all_atoll, num_species, all_nesting, vgb..., all_pc, all_presence);
chain3 = sample(m3, HMC(0.01, 10), 6_000, discard_initial=2000)

using Serialization

serialize("chain_nesting.jls", chain3)

chain3 = deserialize("chain_nesting.jls")

printchain(chain3)

function plotparams(chain::Chains, param::String, n::Int64; lab=nothing)

    # add regex to automatically figure out the n parameters matching param

    p = fill(plot(), n)

    for i in 1:n
        
        sm = @chain DataFrame(group(chain, "$param$i")) begin
            stack(_)
            groupby(_, :variable)
            combine(_, :value => (x -> (mean=mean(logistic.(x)), std=std(logistic.(x)))) => AsTable)
        end

        p[i] = scatter(sm.mean, eachindex(sm.mean), xerror=sm.std, title="$param$i", color=i, label=false)
        
        if !isnothing(lab)
            yt = i in (1, 4) ? lab : fill(" ", length(lab))
            yticks!([0.5:1:length(lab)-0.5;], yt)
            yflip!()
        end

    end

    out = plot([i for i in p]..., size=(1000, 1000))

    return out
end

nestingtypes = String.(unique(Preprocess.envs_known.nestingtype))
species_names = replace.(sort(unique(envs_known.species)), r"[a-z]+_" => " ")

plotparams(chain3, "θ", 6, lab=species_names)
xlims!(0, 1)

# 
# 
# @model function interaction_model(atoll, species, pc, nesting, presence)
#     # Number of groups per predictor
#     a = length(unique(atoll))
#     s = length(unique(species))
#     n = length(unique(nesting))
# 
#     # Priors for atolls
#     ᾱ ~ Normal()
#     τ ~ Exponential(1)
#     α ~ filldist(Normal(ᾱ, τ), a)
#     # Priors for species effect
#     β ~ filldist(Normal(0, 1), s)
#     # Priors for species PC
#     θ1 ~ filldist(Normal(0, 1), s)
#     θ2 ~ filldist(Normal(0, 1), s)
#     θ3 ~ filldist(Normal(0, 1), s)
#     θ4 ~ filldist(Normal(0, 1), s)
#     θ5 ~ filldist(Normal(0, 1), s)
#     θ6 ~ filldist(Normal(0, 1), s)
#     # Priors for nesting PC
#     γ1 ~ filldist(Normal(0, 1), n)
#     γ2 ~ filldist(Normal(0, 1), n)
#     γ3 ~ filldist(Normal(0, 1), n)
#     γ4 ~ filldist(Normal(0, 1), n)
#     γ5 ~ filldist(Normal(0, 1), n)
#     γ6 ~ filldist(Normal(0, 1), n)
# 
# 
#     for i in eachindex(presence)
#         v = logistic(
#             α[atoll[i]] +
#             β[species[i]] +
#             θ1[species[i]] * pc[i, 1] +
#             θ2[species[i]] * pc[i, 2] +
#             θ3[species[i]] * pc[i, 3] +
#             θ4[species[i]] * pc[i, 4] +
#             θ5[species[i]] * pc[i, 5] +
#             θ6[species[i]] * pc[i, 6]
#         )
#         presence[i] ~ Bernoulli(v)
#     end
# end;
# 
# m3 = pc_species_intercept(all_atoll, all_species, all_pc, all_presence)
#chain3 = sample(m3, HMC(0.05, 10), 8000, discard_initial=2000)

# @model function pc_species_intercept(atoll, species, pc, presence)
#     # Number of groups per predictor
#     a = length(unique(atoll))
#     s = length(unique(species))
# 
#     # Priors for atolls
#     # ᾱ ~ Normal()
#     # τ ~ Exponential(1)
#     # α ~ filldist(Normal(ᾱ, τ), a)
#     # Priors for species effect
#     β̄ ~ Normal()
#     κ ~ Exponential(1)
#     β ~ filldist(Normal(β̄, κ), a, s)
#     # Priors for PC
#     θ1 ~ filldist(Normal(0, 1), s)
#     θ2 ~ filldist(Normal(0, 1), s)
#     θ3 ~ filldist(Normal(0, 1), s)
#     θ4 ~ filldist(Normal(0, 1), s)
#     θ5 ~ filldist(Normal(0, 1), s)
#     θ6 ~ filldist(Normal(0, 1), s)
# 
#     for i in eachindex(presence)
#         v = logistic(
#             # α[atoll[i]] +
#             β[atoll[i], species[i]] +
#             θ1[species[i]] * pc[i, 1] +
#             θ2[species[i]] * pc[i, 2] +
#             θ3[species[i]] * pc[i, 3] +
#             θ4[species[i]] * pc[i, 4] +
#             θ5[species[i]] * pc[i, 5] +
#             θ6[species[i]] * pc[i, 6]
#         )
#         presence[i] ~ Bernoulli(v)
#     end
# end;

# m4 = pc_species_intercept(all_atoll, all_species, all_pc, all_presence);
# chain4 = sample(m4, HMC(0.01, 10), 20_000, discard_initial=5000)

# using Serialization

# serialize("chain4.jls", chain4)

# Fit smaller sample
# sdf = @chain Preprocess.envscores begin
#     subset(_, :species => ByRow(x -> x in ("Anous_stolidus", "Phaethon_lepturus", "Onychoprion_fuscatus")), skipmissing=true)
#     # subset(_, :atoll => ByRow(x -> x in ("Tetiaroa", "Palmyra", "Cocos", "Midway", "Wake", "Laysan")), skipmissing=true)
# end
# 
# using StatsBase
# 
# sub_atoll = denserank(sdf.atoll) .|> Int64
# sub_species = denserank(sdf.species) .|> Int64
# sub_pc = Matrix(sdf[:, 1:6])
# sub_presence = sdf.presence .|> Float64
# 
# m4 = pc_species_intercept(sub_atoll, sub_species, sub_pc, sub_presence)
# chain4 = sample(m4, HMC(0.01, 10), 20_000, discard_initial=2000)
# 
# function printchain(chain)
#     x = describe(chain)[1]
#     show(DataFrame(x), allrows=true)
# end
# 
# printchain(chain4)
