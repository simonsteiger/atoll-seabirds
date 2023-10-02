using Turing, StatsFuns, StatsPlots, Chain, DataFrames

include("preprocess.jl")
include("tune.jl")

using .Preprocess
using ReverseDiff

# For saving stuff
using Serialization, Dates

const TODAY = Dates.format(now(), "yyyy-mm-dd")

Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

function printchain(chain)
    x = describe(chain)[1]
    show(DataFrame(x), allrows=true)
end

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

# TODO dont throw out the albatrosses (dont exclude)

@model function pc_species_nesting(atoll, species, nesting, v, g, b, pc, presence)
    # Number of groups per predictor
    a = length(unique(atoll))
    s = length(unique(species))

    # Priors for atolls
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), a)

    # Priors for species effect
    β ~ filldist(Normal(0, 1), s)

    # PC1 per nesting type
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

    # PC2 per nesting type
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

    # PC3 per nesting type
    ω̄31 ~ Normal() # Distribution of distributions of ground nesters
    ω̄32 ~ Normal()
    ω̄33 ~ Normal()
    κ31 ~ Exponential(1)
    κ32 ~ Exponential(1)
    κ33 ~ Exponential(1)
    ω31 ~ filldist(Normal(ω̄31, κ31), b) # Vector of distributions of ground nesters
    ω32 ~ filldist(Normal(ω̄32, κ32), g) # Vector of distributions of tree nesters
    ω33 ~ filldist(Normal(ω̄33, κ33), v) # Vector of distributions of burrow nesters
    ω3 = [ω31, ω32, ω33]

    # PC4 per nesting type
    ω̄41 ~ Normal() # Distribution of distributions of ground nesters
    ω̄42 ~ Normal()
    ω̄43 ~ Normal()
    κ41 ~ Exponential(1)
    κ42 ~ Exponential(1)
    κ43 ~ Exponential(1)
    ω41 ~ filldist(Normal(ω̄41, κ41), b) # Vector of distributions of ground nesters
    ω42 ~ filldist(Normal(ω̄42, κ42), g) # Vector of distributions of tree nesters
    ω43 ~ filldist(Normal(ω̄43, κ43), v) # Vector of distributions of burrow nesters
    ω4 = [ω41, ω42, ω43]
    
    # PC5 per nesting type
    ω̄51 ~ Normal() # Distribution of distributions of ground nesters
    ω̄52 ~ Normal()
    ω̄53 ~ Normal()
    κ51 ~ Exponential(1)
    κ52 ~ Exponential(1)
    κ53 ~ Exponential(1)
    ω51 ~ filldist(Normal(ω̄51, κ51), b) # Vector of distributions of ground nesters
    ω52 ~ filldist(Normal(ω̄52, κ52), g) # Vector of distributions of tree nesters
    ω53 ~ filldist(Normal(ω̄53, κ53), v) # Vector of distributions of burrow nesters
    ω5 = [ω51, ω52, ω53]
    
    # PC6 per nesting type
    ω̄61 ~ Normal() # Distribution of distributions of ground nesters
    ω̄62 ~ Normal()
    ω̄63 ~ Normal()
    κ61 ~ Exponential(1)
    κ62 ~ Exponential(1)
    κ63 ~ Exponential(1)
    ω61 ~ filldist(Normal(ω̄61, κ61), b) # Vector of distributions of ground nesters
    ω62 ~ filldist(Normal(ω̄62, κ62), g) # Vector of distributions of tree nesters
    ω63 ~ filldist(Normal(ω̄63, κ63), v) # Vector of distributions of burrow nesters
    ω6 = [ω61, ω62, ω63]

    for i in eachindex(presence)
        v = logistic(
            α[atoll[i]] +
            β[species[i]] +
            ω1[nesting[i]][species[i]] * pc[i, 1] +
            ω2[nesting[i]][species[i]] * pc[i, 2] +
            ω3[nesting[i]][species[i]] * pc[i, 3] +
            ω4[nesting[i]][species[i]] * pc[i, 4] +
            ω5[nesting[i]][species[i]] * pc[i, 5] +
            ω6[nesting[i]][species[i]] * pc[i, 6]
        )
        presence[i] ~ Bernoulli(v)
    end
end;

m3 = pc_species_nesting(all_atoll, num_species, all_nesting, vgb..., all_pc, all_presence);
chain3 = sample(m3, HMC(0.01, 10), MCMCThreads(), 30_000, 4, discard_initial=5000)

serialize("../chains/$(TODAY)_chain_nesting.jls", chain3)

chain3 = deserialize("chains/chain_nesting.jls")
# 
# printchain(chain3)
# 
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

    out = plot([i for i in p]..., size=(1000, 500), layout=(1,n))

    return out
end

nestingtypes = String.(unique(Preprocess.envs_known.nestingtype))
species_names = replace.(sort(unique(Preprocess.envs_known.species)), r"[a-z]+_" => " ")

plotparams(chain3, "ω6", 3, lab=nothing)
xlims!(0, 1)

df_chain = @chain DataFrame(chain3) begin
    select(_, r"ω̄")
    stack(_)
    groupby(_, :variable)
    combine(_, :value => (x -> (std=std(logistic.(x)), mean=mean(logistic.(x)))) .=> AsTable)
end

[]

groups = chop.(df_chain.variable, tail=1)
scatter(df_chain.mean, df_chain.variable, xerror=df_chain.std, group=groups)
yticks!((1:18).-0.5, df_chain.variable)
xlims!(0, 1)
yflip!()

# Spatial autocorrelation

# using LinearAlgebra
# 
# @model function spatial(atoll, distM, presence)
#     rhosq ~ truncated(Cauchy(0, 1), 0, Inf)
#     etasq ~ truncated(Cauchy(0, 1), 0, Inf)
# 
#     # GPL2
#     SIGMA_distM = etasq * exp.(-rhosq * distM.^2)
#     SIGMA_distM = SIGMA_distM + 0.1I # seems that this step is critical – setting I to 1 seems to avoid the error??
#     SIGMA_distM = (SIGMA_distM' + SIGMA_distM) / 2
#     g ~ MvNormal(zeros(size(SIGMA_distM, 1)), SIGMA_distM)
# 
#     for i in eachindex(presence)
#         lambda = logistic(g[atoll[i]])
#         presence[i] ~ Bernoulli(lambda)
#     end
# end
# 
# m_spatial = spatial(all_atoll, distM_known, all_presence)
# chain_spatial = sample(m_spatial, HMC(0.01, 10), 500)
