using CSV, DataFrames, Chain, Turing, StatsPlots
using MLUtils
using LinearAlgebra

# Where do we get rescale! from?
# Not from MLUtils anyway...
using MLDataUtils: shuffleobs, stratifiedobs, oversample, rescale!

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
    #DataFrames.transform(:nbirds => ByRow(log1p) => :lognbirds)
    subset(_, :nbirds => ByRow(x -> x != 0))
    leftjoin(_, envscores, on=:atoll)
end

# split data function
function split_data(df, species; p=0.50)
    speciesdf = @chain df begin
        subset(_, :species => ByRow(x -> x .== species))
    end
    #shuffled = shuffleobs(speciesdf)
    train, test = splitobs(speciesdf; at=p, shuffle=true) .|> getobs
    # Below code upsamples PC1, but we want to upsample 1-6 simultaneously. Concatenate vectors to array - how?
    # return trainset, testset = oversample(speciesdf.PC1, speciesdf.presence))
end

# Dicts for trainset, testset
trainset = Dict()
testset = Dict()

# Probably fails because some species can't be split without filter conditions
[(trainset[s], testset[s]) = split_data(pop, s) for s in unique(pop.species)];

numerics = Symbol.(PC_NAMES)
features = numerics
target = :nbirds

for f in numerics, k in keys(trainset)
    μ, σ = rescale!(trainset[k][!, f])
    rescale!(testset[k][!, f], μ, σ)
end

# Throw out atolls with mostly 0 counts

# Dicts for train, test, train_label, test_label
train = Dict{String,Matrix}()
test = Dict{String,Matrix}()
train_label = Dict{String,Vector}()
test_label = Dict{String,Vector}()

for k in keys(trainset)
    train[k] = Matrix(trainset[k][:, features])
    test[k] = Matrix(testset[k][:, features])
    train_label[k] = trainset[k][:, target]
    test_label[k] = testset[k][:, target]
end

@model function linear_regression(x, nbirds, n; μ_intercept=0, σ_intercept=1, μ_slope=0, σ_slope=1)
    intercept ~ Normal(μ_intercept, σ_intercept)

    pc1 ~ Normal(μ_slope, σ_slope)
    pc2 ~ Normal(μ_slope, σ_slope)
    pc3 ~ Normal(μ_slope, σ_slope)
    pc4 ~ Normal(μ_slope, σ_slope)
    pc5 ~ Normal(μ_slope, σ_slope)
    pc6 ~ Normal(μ_slope, σ_slope)
    σ ~ Exponential(1)

    for i in 1:n
        μ = intercept + pc1 * x[i, 1] + pc2 * x[i, 2] + pc3 * x[i, 3] + pc4 * x[i, 4] + pc5 * x[i, 5] + pc6 * x[i, 6]
        nbirds[i] ~ Normal(μ, σ)
    end
end;

species = "Anous_minutus"

n, _ = size(train[species])

function standardise(x::AbstractArray)
    x̄, σ = mean(x), std(x)
    [(x - x̄) / σ for x in x]
end

unstandardise(x::AbstractArray, x̄::Real, σ::Real) = x * σ .+ x̄

m = linear_regression(train[species], standardise(log.(train_label[species])), n) # log train label
chain = sample(m, NUTS(), MCMCThreads(), 30_000, 3)

params = select(DataFrame(chain), r"intercept|pc")

plot(chain)

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
end

out = prediction(test[species], chain)

x̄ = mean(log.(test_label[species]))
σ = std(log.(test_label[species]))

natural = unstandardise(out, x̄, σ)

Q = [quantile(slice, limits) for limits in [0.025, 0.975], slice in eachslice(natural, dims=1)]

μ_natural = [mean(x) for x in eachslice(natural, dims=1)]

mean(exp.(μ_natural))
std(exp.(μ_natural))
mean(test_label[species])
std(test_label[species])

sum(exp.(μ_natural))
sum(test_label[species])

scatter(μ_natural, label="predicted", yerror=Q, alpha=0.5)
scatter!(log.(test_label[species]), label="observed", alpha=0.5)
#ylims!((-10, 1e4))
