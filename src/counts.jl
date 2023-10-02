module Counts

using CSV, DataFrames, Chain, Turing, StatsPlots
using MLUtils
using LinearAlgebra

# Where do we get rescale! from?
# Not from MLUtils anyway...
using MLDataUtils: shuffleobs, stratifiedobs, oversample, rescale!

include("standardise.jl")

const PC_NAMES = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"]

envscores = @chain begin
    CSV.read("data/jl_envscores.csv", DataFrame)
    select(_, [:atoll], Symbol.(PC_NAMES))
    unique(_, [:atoll])
end

pop = @chain begin
    CSV.read("data/atoll_seabird_populations_29Jul.csv", DataFrame)
    DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    subset(_, All() .=> ByRow(x -> x != -1))
    stack(_, Not(:atoll), variable_name=:species, value_name=:nbirds)
    subset(_, :nbirds => ByRow(x -> x != 0))
    leftjoin(_, envscores, on=:atoll)
end

numerics = Symbol.(PC_NAMES)
features = numerics
target = :nbirds

spdict = Dict()

for species in unique(pop.species)
    spdict[species] = @chain pop begin
        subset(_, :species => ByRow(x -> x .== species))
    end
end

spfeatures = Dict()
sptarget = Dict()

for k in keys(spdict)
    spfeatures[k] = Matrix(spdict[k][:, features])
    sptarget[k] = spdict[k][:, target]
end

@model function linear_regression(x, nbirds, n; df=3)
    intercept ~ TDist(df)

    pc1 ~ TDist(df)
    pc2 ~ TDist(df)
    pc3 ~ TDist(df)
    pc4 ~ TDist(df)
    pc5 ~ TDist(df)
    pc6 ~ TDist(df)
    σ ~ Exponential(1)

    for i in 1:n
        μ = intercept + pc1 * x[i, 1] + pc2 * x[i, 2] + pc3 * x[i, 3] + pc4 * x[i, 4] + pc5 * x[i, 5] + pc6 * x[i, 6]
        nbirds[i] ~ Normal(μ, σ)
    end
end;

function prediction(x, chain)
    n, _ = size(x)

    params = select(DataFrame(chain), r"pc|intercept")
    chainlength = nrow(params)
    intercept = params.intercept
    pc1 = params.pc1
    pc2 = params.pc2
    pc3 = params.pc3
    pc4 = params.pc4
    pc5 = params.pc5
    pc6 = params.pc6

    nbirds = zeros(n, chainlength)

    for i in 1:n, j in 1:chainlength
        nbirds[i, j] = intercept[j] + pc1[j] * x[i, 1] + pc2[j] * x[i, 2] + pc3[j] * x[i, 3] + pc4[j] * x[i, 4] + pc5[j] * x[i, 5] + pc6[j] * x[i, 6]
    end

    return nbirds
end;

struct result
    chain::Chains
    df::Matrix
end

function wrap_model(features::Dict, target::Dict, species)
    n, _ = size(features[species])
    m = linear_regression(features[species], standardise(log10.(target[species])), n)
    chain = sample(m, NUTS(), MCMCThreads(), 30_000, 3)
    df = prediction(features[species], chain)
    result(chain, df)
end;

# FIX 
# Take away species with insufficient data
# [wrap_model(spfeatures, sptarget, s) for s in keys(spfeatures)]

species = "Anous_stolidus"

out = wrap_model(spfeatures, sptarget, species)

x̄ = mean(log10.(sptarget[species]))
σ = std(log10.(sptarget[species]))

natural = unstandardise(out.df, x̄, σ)

Q = [quantile(row, limits) for limits in [0.005, 0.995], row in eachslice(natural, dims=1)]

μ_natural = [mean(x) for x in eachslice(natural, dims=1)]

scatter(log10.(sptarget[species]), label="observed", color="black")
scatter!(μ_natural, label="predicted", marker=:d, color="red", yerror=(abs.(μ_natural.-Q[1, :]), abs.(μ_natural.-Q[2, :])))
title!(replace(species, Pair("_", " ")))
ylims!((0, Inf))

end

species = "Anous_stolidus"

n, _ = size(spfeatures[species])

m = linear_regression(spfeatures[species], standardise(log10.(sptarget[species])), n); # log train label
chain = sample(m, NUTS(), MCMCThreads(), 30_000, 3)

params = select(DataFrame(chain), r"intercept|pc")

plot(chain)


out = prediction(spfeatures[species], chain)

x̄ = mean(log10.(sptarget[species]))
σ = std(log10.(sptarget[species]))

natural = unstandardise(out, x̄, σ)

Q = [quantile(row, limits) for limits in [0.005, 0.995], row in eachslice(natural, dims=1)]

μ_natural = [mean(x) for x in eachslice(natural, dims=1)]

mean(exp10.(μ_natural))
std(exp10.(μ_natural))
mean(sptarget[species])
std(sptarget[species])

sum(exp10.(μ_natural))
sum(sptarget[species])

scatter(μ_natural, label="predicted", yerror=transpose(Q), alpha=0.5)
scatter!(log10.(sptarget[species]), label="observed", alpha=0.5)
title!(species)
ylims!((0, Inf))
