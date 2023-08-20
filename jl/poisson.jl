using CSV, DataFrames, Chain, Turing, StatsPlots
using MLUtils

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

let v = pop.PC1
    μ = mean(v)
    σ = std(v)
    histogram((v .- μ) / σ)
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

@model function poisson_regression(x, nbirds, n; μ_intercept=3, σ_intercept=0.5, μ_slope=0, σ_slope=0.1)
    intercept ~ Normal(μ_intercept, σ_intercept)

    pc1 ~ Normal(μ_slope, σ_slope)
    pc2 ~ Normal(μ_slope, σ_slope)
    pc3 ~ Normal(μ_slope, σ_slope)
    pc4 ~ Normal(μ_slope, σ_slope)
    pc5 ~ Normal(μ_slope, σ_slope)
    pc6 ~ Normal(μ_slope, σ_slope)
    p ~ Beta(2, 2)

    for i in 1:n
        λ = exp(intercept + pc1 * x[i, 1] + pc2 * x[i, 2] + pc3 * x[i, 3] + pc4 * x[i, 4] + pc5 * x[i, 5] + pc6 * x[i, 6])
        nbirds[i] ~ NegativeBinomial(λ+1e-5, p)
    end
end;

species = "Gygis_alba"

n, _ = size(train[species])

m = poisson_regression(train[species], train_label[species], n)
chain = sample(m, NUTS(), MCMCThreads(), 10_000, 3)

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
        nbirds[i, j] = exp(intercept[j] + pc1[j] * x[i, 1] + pc2[j] * x[i, 2] + pc3[j] * x[i, 3] + pc4[j] * x[i, 4] + pc5[j] * x[i, 5] + pc6[j] * x[i, 6])
    end

    return nbirds
end

out = prediction(test[species], chain)

out_mean = [mean(x) for x in eachslice(out, dims=1)]

mean(out_mean)
mean(test_label[species])

scatter(out_mean, label = "predicted")
scatter!(test_label[species], label = "observed")
#ylims!((-10, 1e4))
