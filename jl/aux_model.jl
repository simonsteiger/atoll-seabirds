include("/Users/simonsteiger/Desktop/other/atoll-seabirds/jl/aux_prep.jl")

using Turing

# Bayesian logistic regression with Student T priors
@model function logistic_regression(x, y, n, df)
    intercept ~ TDist(df)

    pc1 ~ TDist(df)
    pc2 ~ TDist(df)
    pc3 ~ TDist(df)
    pc4 ~ TDist(df)
    pc5 ~ TDist(df)
    pc6 ~ TDist(df)

    for i in 1:n
        v = logistic(intercept + pc1 * x[i, 1] + pc2 * x[i, 2] + pc3 * x[i, 3] + pc4 * x[i, 4] + pc5 * x[i, 5] + pc6 * x[i, 6])
        y[i] ~ Bernoulli(v)
    end
end;


function wrap_model(species, fun, σ; print=false)
    # Retrieve the number of observations
    n, _ = size(train[species])

    # Sample using HMC
    m = fun(train[species], train_label[species], n, σ)
    chain = sample(m, HMC(0.05, 10), MCMCThreads(), 10_000, 3, discard_initial=5000)

    tuning_params = tune(chain, test, test_label)

    tuning_params_sub = @chain tuning_params begin
        subset(_, :species => x -> x .== species)
    end

    optimal = tuning_params_sub.threshold[maximum(tuning_params_sub.criterion).==tuning_params_sub.criterion][1]

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

    diagnostics = Dict(
        "present" => present,
        "absent" => absent,
        "predicted_present" => predicted_present,
        "predicted_absent" => predicted_absent
    )

    if print
        println("
        Species: $species
        ----
        Predicted absent: $predicted_absent
        Observed absent: $absent
        Percentage absent correct $(predicted_absent/absent)
        ----
        Predicted present: $predicted_present
        Observed present: $present
        Percentage present correct $(predicted_present/present)
        ")
    end

    return result(chain, threshold, diagnostics)
end

# Predict binary outcomes from a logistic model, using the mean of each paramter
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
end

# Predict binary outcomes from a logistic regression, using all the sames from each chain
function prediction_full(x::Matrix, chain, threshold)
    # Pull the means from each parameter's sampled values in the chain
    intercept = chain[:intercept]
    pc1 = chain[:pc1]
    pc2 = chain[:pc2]
    pc3 = chain[:pc3]
    pc4 = chain[:pc4]
    pc5 = chain[:pc5]
    pc6 = chain[:pc6]

    # Retrieve number of rows
    n, _ = size(x)

    # Generate a vector to store our predictions
    v = Matrix{Float64}(undef, n, length(intercept))

    # Calculate the logistic function for each element in the test set
    for i in 1:n, j in 1:30_000
        num = logistic(
            intercept[j] + pc1[j] * x[i, 1] + pc2[j] * x[i, 2] + pc3[j] * x[i, 3] + pc4[j] * x[i, 4] + pc5[j] * x[i, 5] + pc6[j] * x[i, 6]
        )
        if num .>= threshold
            v[i, j] = 1
        else
            v[i, j] = 0
        end
    end

    return v
end
