# Import Turing and Distributions.
using Turing, Distributions

# Working with tabular data
using CSV, DataFrames, Chain

# Import MCMCChains, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots

# We need a logistic function, which is provided by StatsFuns.
using StatsFuns: logistic

using StatsBase

# Functionality for splitting and normalizing the data
using MLDataUtils: shuffleobs, stratifiedobs, oversample, rescale!

# Set a seed for reproducibility.
using Random

include("/Users/simonsteiger/Desktop/other/atoll-seabirds/R/upsample.jl")
include("/Users/simonsteiger/Desktop/other/atoll-seabirds/R/tune.jl")

Random.seed!(0);

# Turn off progress monitor.
# Turing.setprogress!(false)

# Import the data
envscores = CSV.read("data/envscores.csv", DataFrame)

# Show first few rows
first(envscores, 5)

# Convert "species" and "presence" to numeric values
envscores[!, :presence] = [r.presence == true ? 1.0 : 0.0 for r in eachrow(envscores)]

# Check if conversion worked
first(envscores, 5)

# Delete all species with known population from data
envscores = subset(envscores, [:cond, :region] => ByRow((x, y) -> occursin(Regex(x), y)))

# Discard unused columns
select!(envscores, Not([:Column1, :region, :cond]))

# split data function
function split_data(df, target, species; at=0.70)
    speciesdf = @chain df begin 
        subset(_, :species => ByRow(x -> x .== species))
    end
    shuffled = shuffleobs(speciesdf)
    return trainset, testset = stratifiedobs(row -> row[target], shuffled; p = at)
    # Below code upsamples PC1, but we want to upsample 1-6 simultaneously. Concatenate vectors to array - how?
    # return trainset, testset = oversample(speciesdf.PC1, speciesdf.presence))
end

features = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
numerics = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
target = :presence

# Dicts for trainset, testset
trainset = Dict{String, DataFrame}()
testset = Dict{String, DataFrame}()

[(trainset[s], testset[s]) = split_data(envscores, target, s; at=0.5) for s in unique(envscores.species)]

for f in numerics, k in keys(trainset)
    μ, σ = rescale!(trainset[k][!, f]; obsdim=1)
    rescale!(testset[k][!, f], μ, σ; obsdim=1)
end

trainset_up = Dict{String, DataFrame}()

[trainset_up[k] = upsample(trainset[k], 250) for k in keys(trainset)];

# Dicts for train, test, train_label, test_label
train = Dict{String,Matrix}()
test = Dict{String,Matrix}()
train_label = Dict{String,Vector}()
test_label = Dict{String,Vector}()

for k in keys(trainset)
    train[k] = Matrix(trainset_up[k][:, features])
    test[k] = Matrix(testset[k][:, features])
    train_label[k] = trainset_up[k][:, target]
    test_label[k] = testset[k][:, target];
end

# Bayesian logistic regression
@model function logistic_regression(x, y, n, σ)
    intercept ~ Normal(0, σ)

    pc1 ~ Normal(0, σ)
    pc2 ~ Normal(0, σ)
    pc3 ~ Normal(0, σ)
    pc4 ~ Normal(0, σ)
    pc5 ~ Normal(0, σ)
    pc6 ~ Normal(0, σ)

    for i in 1:n
        v = logistic(intercept + pc1 * x[i, 1] + pc2 * x[i, 2] + pc3 * x[i, 3] + pc4 * x[i, 4] + pc5 * x[i, 5] + pc6 * x[i, 6])
        y[i] ~ Bernoulli(v)
    end
end;

# As long as pipeline is only partially looped, specify species name here
species = "Anous_stolidus"

# Retrieve the number of observations
n, _ = size(train[species])

# Sample using HMC
m = logistic_regression(train[species], train_label[species], n, 1)
chain = sample(m, HMC(0.05, 10), MCMCThreads(), 10_000, 3, discard_initial = 5000)

# Plot
plot(chain)

# Labels
l = [:pc1, :pc2, :pc3, :pc4, :pc5, :pc6]

# Corner plot
corner(chain, l)

function prediction(x::Matrix, chain, threshold)
    # Pull the means from each parameter's sampled values in the chain
    intercept = mean(chain[:intercept])
    pc1 = mean(chain[:pc1])
    pc2 = mean(chain[:pc2])
    pc3 = mean(chain[:pc3])
    pc4 = mean(chain[:pc4])
    pc5 = mean(chain[:pc5])
    pc6 = mean(chain[:pc6])

    # Retrieve number of rows
    n, _ = size(x)

    # Generate a vector to store out predictions
    v = Vector{Float64}(undef, n)

    # Calculate the logistic function for each element in the test set
    for i in 1:n
        num = logistic(
            intercept .+ pc1 * x[i, 1] + pc2 * x[i, 2] + pc3 * x[i, 3] + pc4 * x[i, 4] + pc5 * x[i, 5] + pc6 * x[i, 6]
        )
        if num >= threshold
            v[i] = 1
        else
            v[i] = 0
        end
    end
    return v
end;

tuning_params = tune(chain, test, test_label)

tuning_params_sub = @chain tuning_params begin
    subset(_, :species => x -> x .== species)
end

plot(tuning_params_sub.threshold, tuning_params_sub.criterion)
plot!(tuning_params_sub.threshold, tuning_params_sub.predicted_absent)
plot!(tuning_params_sub.threshold, tuning_params_sub.predicted_present)

optimal = tuning_params_sub.threshold[maximum(tuning_params_sub.criterion) .== tuning_params_sub.criterion][1]

scatter!([optimal], [maximum(tuning_params_sub.criterion)])

# Set predictions threshold
threshold = optimal

# Make the predictions
predictions = prediction(test[species], chain, threshold)

# Calculate MSE for our test set
loss = sum((predictions - test_label[species]) .^ 2) / length(test_label[species])

present = sum(test_label[species])
absent = length(test_label[species]) - present

predicted_present = sum(test_label[species] .== predictions .== 1)
predicted_absent = sum(test_label[species] .== predictions .== 0)

println("Predicted absent: $predicted_absent
Observed absent: $absent
Percentage absent correct $(predicted_absent/absent)
---
Predicted present: $predicted_present
Observed present: $present
Percentage present correct $(predicted_present/present)")
