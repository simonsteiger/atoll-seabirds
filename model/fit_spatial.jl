using Turing, StatsFuns, StatsPlots, Chain, DataFrames

include("../src/preprocess.jl")
include("../src/tune.jl")

using .Preprocess
using ReverseDiff
using StatsBase
using LinearAlgebra

# For saving stuff
using Serialization, Dates

const TODAY = Dates.format(now(), "yyyy-mm-dd")

# Should the model be fit?
fit = false

Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

function printchain(chain)
    x = describe(chain)[1]
    show(DataFrame(x), allrows=true)
end


@model function full_model(atoll, species, nesting, v, g, b, pc, dist, presence)
    # Number of groups per predictor
    a = length(unique(atoll))
    s = length(unique(species))

    # Priors for atolls
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), a)

    # Priors for species effect
    β ~ filldist(Normal(0, 1), s)

    # Spatial stuff
    rhosq ~ truncated(Cauchy(0, 1), 0, Inf)
    etasq ~ truncated(Cauchy(0, 1), 0, Inf)

    # GPL2
    SIGMA_distM = etasq * exp.(-rhosq * distM.^2)
    SIGMA_distM = SIGMA_distM + 0.01I # seems that this step is critical – setting I to 1 seems to avoid the error??
    SIGMA_distM = (SIGMA_distM' + SIGMA_distM) / 2
    g ~ MvNormal(zeros(size(SIGMA_distM, 1)), SIGMA_distM)

    for i in eachindex(presence)
	    lambda = logistic(α[atoll[i]] + β[species[i]] +g[atoll[i]])
        presence[i] ~ Bernoulli(lambda)
    end
end

m_spatial = spatial(all_atoll, distM_known, all_presence)

if fit
    chain_spatial = sample(m_spatial, HMC(0.01, 10), 5000, discard_initial=1000)
    serialize("spatial_chain.jls", chain_spatial)
end

prior_samples = sample(m_spatial, Prior(), 1000)

df_ps = DataFrame(prior_samples)
select!(df_ps, [:rhosq, :etasq])

p = plot()

f(ρ, η) = [η * exp(-ρ * x^2) for x in 0:0.1:15]
plot!(f.(df_ps.rhosq, df_ps.etasq), legend=false)
xticks!([0:10:40;], string.([0:4;]))
xlabel!("distance (1000 km)")

using Serialization

chain_spatial = deserialize("model/chains/2023-10-07_spatial_chain.jls")

df_cs = DataFrame(chain_spatial)
select!(df_cs, [:rhosq, :etasq])

plot!(f.(df_cs.rhosq[begin:40:end], df_cs.etasq[begin:40:end]), color=2, legend=false, alpha=0.2)
xticks!([0:10:150;], string.([0:15;]))
xlabel!("distance (1000 km)")

