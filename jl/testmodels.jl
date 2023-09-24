using Turing, StatsFuns, StatsPlots, Chain, DataFrames

include("preprocess.jl")
include("tune.jl")

using .Preprocess, ReverseDiff

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

@model function pc_species_nesting(atoll, species, nesting, pc, presence)
    # Number of groups per predictor
    a = length(unique(atoll))
    s = length(unique(species))
    n = length(unique(nesting))

    # Priors for atolls
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), a)

    # Priors for species effect
    β ~ filldist(Normal(0, 1), s)
    
    # Priors for species-PC
    θ1 ~ filldist(Normal(0, 1), s)
    θ2 ~ filldist(Normal(0, 1), s)
    θ3 ~ filldist(Normal(0, 1), s)
    θ4 ~ filldist(Normal(0, 1), s)
    θ5 ~ filldist(Normal(0, 1), s)
    θ6 ~ filldist(Normal(0, 1), s)
    
    # Priors for nesting-PC
    γ̄1 ~ Normal()
    κ1 ~ Exponential(1)
    γ1 ~ filldist(Normal(γ̄1, κ1), n)
    
    γ̄2 ~ Normal()
    κ2 ~ Exponential(1)
    γ2 ~ filldist(Normal(γ̄2, κ2), n)
    
    γ̄3 ~ Normal()
    κ3 ~ Exponential(1)
    γ3 ~ filldist(Normal(γ̄3, κ3), n)
    
    γ̄4 ~ Normal()
    κ4 ~ Exponential(1)
    γ4 ~ filldist(Normal(γ̄4, κ4), n)
    
    γ̄5 ~ Normal()
    κ5 ~ Exponential(1)
    γ5 ~ filldist(Normal(γ̄5, κ5), n)
    
    γ̄6 ~ Normal()
    κ6 ~ Exponential(1)
    γ6 ~ filldist(Normal(γ̄6, κ6), n)

    for i in eachindex(presence)
        v = logistic(
            α[atoll[i]] +
            β[species[i]] +
            θ1[species[i]] * pc[i, 1] +
            θ2[species[i]] * pc[i, 2] +
            θ3[species[i]] * pc[i, 3] +
            θ4[species[i]] * pc[i, 4] +
            θ5[species[i]] * pc[i, 5] +
            θ6[species[i]] * pc[i, 6] +
            γ1[nesting[i]] * pc[i, 1] +
            γ2[nesting[i]] * pc[i, 2] +
            γ3[nesting[i]] * pc[i, 3] +
            γ4[nesting[i]] * pc[i, 4] +
            γ5[nesting[i]] * pc[i, 5] +
            γ6[nesting[i]] * pc[i, 6]
        )
        presence[i] ~ Bernoulli(v)
    end
end;

m3 = pc_species_nesting(all_atoll, all_species, all_nesting, all_pc, all_presence);
chain3 = sample(m3, HMC(0.01, 10), 15_000, discard_initial=5000)

using Serialization

serialize("chain_nesting.jls", chain3)

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

plotparams(chain3, "γ", 6, nestingtypes)
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
