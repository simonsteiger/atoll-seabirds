module MultiPresence

using Turing, StatsFuns, StatsPlots

include("preprocess.jl")
include("tune.jl")

using .Preprocess

# Multilevel model
# Hmmm ... shouldn't we have a species intercept?
@model function multilevelmodel(pc, presence, atoll)
    k = length(unique(atoll))
    ᾱ ~ Normal()
    τ ~ Exponential(1)
    α ~ filldist(Normal(ᾱ, τ), k)

    β ~ filldist(Normal(0, 1), 6)

    for i in eachindex(presence)
        v = logistic(α[atoll[i]] + β[1] * pc[i, 1] + β[2] * pc[i, 2] + β[3] * pc[i, 3] + β[4] * pc[i, 4] + β[5] * pc[i, 5] + β[6] * pc[i, 6])
        presence[i] ~ Bernoulli(v)
    end
end;

s = "Anous_stolidus"
m = multilevelmodel(train[s], train_label[s], groupatoll[s])
chain = sample(m, HMC(0.05, 10), MCMCThreads(), 15_000, 4, discard_initial=5000)

plot(logistic.(mean(group(chain, "α"))[:, 2]))

function prediction(x::Matrix, chain, threshold)
    # Retrieve number of rows
    n, _ = size(x)

    # Generate a vector to store out predictions
    v = Vector{Float64}(undef, n)

    # Calculate the logistic function for each element in the test set
    for i in 1:n
        num = logistic(
            mean(chain["α[$i]"]) .+ mean(chain["β[1]"]) * x[i, 1]
            + mean(chain["β[2]"]) * x[i, 2]
            + mean(chain["β[3]"]) * x[i, 3]
            + mean(chain["β[4]"]) * x[i, 4]
            + mean(chain["β[5]"]) * x[i, 5]
            + mean(chain["β[6]"]) * x[i, 6]
        )
        if num >= threshold
            v[i] = 1
        else
            v[i] = 0
        end
    end
    return v
end

tuning_params = tune(chain, Preprocess.test[s], Preprocess.test_label[s])

optimal = tuning_params.threshold[maximum(tuning_params.criterion).==tuning_params.criterion][1]

# scatter!([optimal], [maximum(tuning_params_sub.criterion)])

# Set predictions threshold
threshold = optimal

# Make the predictions
predictions = prediction(Preprocess.test[s], chain, threshold)

present = sum(test_label[s])
absent = length(test_label[s]) - present

predicted_present = sum(test_label[s] .== predictions .== 1)
predicted_absent = sum(test_label[s] .== predictions .== 0)

predicted_present / present
predicted_absent / absent

end