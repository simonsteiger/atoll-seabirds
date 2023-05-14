# Import Turing and Distributions.
using Turing, Distributions

# Import CSV
using CSV

# Import MCMCChains, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots

# We need a logistic function, which is provided by StatsFuns.
using StatsFuns: logistic

# Functionality for splitting and normalizing the data
using MLDataUtils: shuffleobs, stratifiedobs, rescale!

# Set a seed for reproducibility.
using Random
Random.seed!(0);

# Turn off progress monitor.
# Turing.setprogress!(false)

# Import the data
envscores = CSV.read("data/envscores.csv", DataFrame)

# Show first few rows
first(envscores, 5)

# Convert "species" and "presence" to numeric values
envscores[!, :presence] = [r.presence == true ? 1.0 : 0.0 for r in eachrow(envscores)]
# envscores[!, :species] = indexin(envscores[!, :species], unique(envscores.species))
# envscores[!, :species] = Float64.(envscores[!, :species])
select!(envscores, Not(:Column1))
subset!(envscores, :species => ByRow(x -> x .== "Hydroprogne_caspia"))


# Check if conversion worked
first(envscores, 5)

# split data function
function split_data(df, target; at=0.70)
    shuffled = shuffleobs(df)
    return trainset, testset = stratifiedobs(row -> row[target], shuffled; p = at)
end

features = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
numerics = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
target = :presence

trainset, testset = split_data(envscores, target; at=0.5)

for f in numerics
    μ, σ = rescale!(trainset[!, f]; obsdim=1)
    rescale!(testset[!, f], μ, σ; obsdim=1)
end

train = Matrix(trainset[:, features])
test = Matrix(testset[:, features])
train_label = trainset[:, target]
test_label = testset[:, target];

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

# Retrieve the number of observations
n, _ = size(train)

# Sample using HMC
m = logistic_regression(train, train_label, n, 1)
chain = sample(m, HMC(0.05, 10), MCMCThreads(), 10_000, 3)

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

# Set predictions threshold
threshold = 0.225

# Make the predictions
predictions = prediction(test, chain, threshold)

# Calculate MSE for our test set
loss = sum((predictions - test_label) .^ 2) / length(test_label)

present = sum(test_label)
absent = length(test_label) - present

predicted_present = sum(test_label .== predictions .== 1)
predicted_absent = sum(test_label .== predictions .== 0)

println("Predicted absent: $predicted_absent
Observed absent: $absent
Percentage absent correct $(predicted_absent/absent)
---
Predicted present: $predicted_present
Observed present: $present
Percentage present correct $(predicted_present/present)")
