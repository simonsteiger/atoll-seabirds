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
using Resample

# Set a seed for reproducibility.
using Random, StableRNGs

#import .Upsample

include("/Users/simonsteiger/Desktop/other/atoll-seabirds/jl/upsample.jl")
include("/Users/simonsteiger/Desktop/other/atoll-seabirds/jl/tune.jl")

Random.seed!(0);

const PRIOR_TEST = 0.5
const PRIOR_DF = 3

# Turn off progress monitor.
# Turing.setprogress!(false)

# Import the data
envscores = CSV.read("data/jl_envscores.csv", DataFrame)

# Convert "presence" to numeric values
# envscores[!, :presence] = [r.presence == true ? 1.0 : 0.0 for r in eachrow(envscores)]

envs_known = subset(envscores, :presence => ByRow(x -> !ismissing(x)))
envs_unknown = subset(envscores, :presence => ByRow(x -> ismissing(x)))

# Delete all species with known population from data
envs_known = subset(envs_known, [:filtercondition, :region] => ByRow((x, y) -> occursin(Regex(x), y)))

# Discard unused columns
select!(envs_known, Not([:region, :filtercondition]))

# split data function
function split_data(df, target, species; at=0.70)
    speciesdf = @chain df begin
        subset(_, :species => ByRow(x -> x .== species))
    end
    shuffled = shuffleobs(speciesdf)
    return trainset, testset = stratifiedobs(row -> row[target], shuffled; p=at)
    # Below code upsamples PC1, but we want to upsample 1-6 simultaneously. Concatenate vectors to array - how?
    # return trainset, testset = oversample(speciesdf.PC1, speciesdf.presence))
end

features = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
numerics = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]
target = :presence

# Dicts for trainset, testset
trainset = Dict{String,DataFrame}()
testset = Dict{String,DataFrame}()

[(trainset[s], testset[s]) = split_data(envs_known, target, s; at=0.5) for s in unique(envs_known.species)]

for f in numerics, k in keys(trainset)
    μ, σ = rescale!(trainset[k][!, f]; obsdim=1)
    rescale!(testset[k][!, f], μ, σ; obsdim=1)
end

trainset_up = Dict{String,DataFrame}()

rng = StableRNG(1)

[trainset_up[k] = our_smote(rng, trainset[k]) for k in keys(trainset)];

# Dicts for train, test, train_label, test_label
train = Dict{String,Matrix}()
test = Dict{String,Matrix}()
train_label = Dict{String,Vector}()
test_label = Dict{String,Vector}()

for k in keys(trainset)
    train[k] = Matrix(trainset_up[k][:, features])
    test[k] = Matrix(testset[k][:, features])
    train_label[k] = trainset_up[k][:, target]
    test_label[k] = testset[k][:, target]
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

# Bayesian logistic regression
@model function logistic_regressionT(x, y, n, df)
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

struct result
    chain::Chains
    threshold::Float64
    diagnostics::Dict{String,Integer}
end

function wrap_model(species, model_fun, dispersion)
    # Retrieve the number of observations
    n, _ = size(train[species])

    # Sample using HMC
    m = model_fun(train[species], train_label[species], n, dispersion)
    chain = sample(m, HMC(0.05, 10), MCMCThreads(), 10_000, 3, discard_initial=5000)

    # Plot
    plot(chain)

    # Labels
    l = [:pc1, :pc2, :pc3, :pc4, :pc5, :pc6]

    # Corner plot
    # corner(chain, l)

    tuning_params = tune(chain, test, test_label)

    tuning_params_sub = @chain tuning_params begin
        subset(_, :species => x -> x .== species)
    end

    # plot(tuning_params_sub.threshold, tuning_params_sub.criterion)
    # plot!(tuning_params_sub.threshold, tuning_params_sub.predicted_absent)
    # plot!(tuning_params_sub.threshold, tuning_params_sub.predicted_present)

    optimal = tuning_params_sub.threshold[maximum(tuning_params_sub.criterion).==tuning_params_sub.criterion][1]

    # scatter!([optimal], [maximum(tuning_params_sub.criterion)])

    # Set predictions threshold
    threshold = optimal

    # Make the predictions
    predictions = prediction(test[species], chain, threshold)

    # Calculate MSE for our test set
    # loss = sum((predictions - test_label[species]) .^ 2) / length(test_label[species])

    # Creating this dictionary could be a function
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

    return result(chain, threshold, diagnostics)
end

# all_results_N = Dict{String,result}()
all_results_T = Dict{String,result}()

[all_results_T[k] = wrap_model(k, logistic_regressionT, PRIOR_DF) for k in collect(keys(trainset))]

# [all_results_N[k] = wrap_model(k, logistic_regression, p) for k in collect(keys(trainset)), p in PRIOR_TEST]

missing_atolls = @chain envs_known begin
    subset(_, :presence => ByRow(x -> ismissing(x)))
end

some_dict_T = Dict()
# some_dict_N = Dict()
X_envs_unknown = Matrix(envs_unknown[!, begin:6])

# [some_dict_N[k] = prediction(X_envs_unknown, all_results_N[k].chain, all_results_N[k].threshold) for k in keys(all_results_N)]

# temp = map([some_dict_T, some_dict_N]) do data
#     @chain DataFrame(data) begin
#         insertcols(_, :atoll => envs_unknown.atoll)
#         stack(_, Not(:atoll), variable_name=:species, value_name=:presence)
#     end
# end

temp = @chain DataFrame(some_dict_T) begin
    insertcols(_, :atoll => envs_unknown.atoll)
    stack(_, Not(:atoll), variable_name=:species, value_name=:presence)
end

# CSV.write("data/first_preds.csv", temp)

# TODO
# Check if there is a correlation between low majority class N and prediction performance
# Can we fit a single model with indices for each species?

function squash(x::AbstractArray)
    return reduce(hcat, x)'
end

function prediction_full(x::Matrix, chain, threshold)
    # Pull the means from each parameter's sampled values in the chain
    intercept = squash(chain[:intercept])
    pc1 = squash(chain[:pc1])
    pc2 = squash(chain[:pc2])
    pc3 = squash(chain[:pc3])
    pc4 = squash(chain[:pc4])
    pc5 = squash(chain[:pc5])
    pc6 = squash(chain[:pc6])

    # Retrieve number of rows
    n, _ = size(x)

    # Generate a vector to store out predictions
    v = Matrix{Float64}(undef, n, length(intercept))
    # num = Vector{Float64}(undef, n)

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

[some_dict_T[k] = prediction_full(X_envs_unknown, all_results_T[k].chain, all_results_T[k].threshold) for k in keys(all_results_T)]

pred_conf = Dict()

map(collect(keys(some_dict_T))) do k
    pred_conf[k] =
        @chain DataFrame(some_dict_T[k], :auto) begin
            insertcols(_, 1, :atoll => Vector{String}(envs_unknown.atoll))
            stack(_)
            groupby(_, :atoll)
            combine(_, :value .=> (x -> mean(x)) => :percent_present)
        end
end

function mini_summary(k; cutoff=0.9)
    pct = sum(pred_conf[k].percent_present .> cutoff) / nrow(pred_conf[k])
    println("""
    $k: $(round(100*pct, digits=1))
    """
    )
end

[mini_summary(k, cutoff=0.8) for k in keys(pred_conf)]

N_pred, T_pred = Dict(), Dict()

# [N_pred[k] = prediction_full(test[k], all_results_N[k].chain, all_results_N[k].threshold)' for k in keys(all_results_N)]

[T_pred[k] = prediction_full(test[k], all_results_T[k].chain, all_results_T[k].threshold)' for k in keys(all_results_T)]

pred_dct = Dict()

# for i in ["N", "T"], k in keys(all_results_N)
#     data = i == "N" ? N_pred : T_pred
#     pred_dct["$k,$i"] =
#         @chain DataFrame(data[k], :auto) begin
#             combine(_, Cols(startswith("x")) .=> (x -> abs(mean(x) - 0.5) * 200))
#             stack(_)
#             # DataFrames.transform(_, :value => ByRow(x -> ≈(x, 100, atol=5)) => :conf)
#         end
# end

p_dict = Dict()

for k in keys(all_results_N)
    p = density()
    [density!(pred_dct["$k,$i"].value, lw=1.5, label=i) for i in ["N", "T"]]
    title!(k)
    p_dict[k] = p
end

# TODO run this plot for wide, default, and narrow settings – e.g. Normal(0, 0.5), TDist(3), TDist(1)
plot([p_dict[k] for k in keys(p_dict)]..., size=(1600, 800), titlefontsize=10)

# Prior predictive checks
# Looks like TDist(3) is a good regularising prior
@model function sample_priorT(x, y, df)
    intercept ~ TDist(df)

    pc1 ~ TDist(df)
    pc2 ~ TDist(df)
    pc3 ~ TDist(df)
    pc4 ~ TDist(df)
    pc5 ~ TDist(df)
    pc6 ~ TDist(df)

    v = logistic(intercept + pc1 + pc2 + pc3 + pc4 + pc5 + pc6)
    y ~ Bernoulli(v)
end;

@model function sample_priorN(x, y, σ)
    intercept ~ Normal(0, σ)

    pc1 ~ Normal(0, σ)
    pc2 ~ Normal(0, σ)
    pc3 ~ Normal(0, σ)
    pc4 ~ Normal(0, σ)
    pc5 ~ Normal(0, σ)
    pc6 ~ Normal(0, σ)

    v = logistic(intercept + pc1 + pc2 + pc3 + pc4 + pc5 + pc6)
    y ~ Bernoulli(v)
end;

pT = sample(sample_priorT(missing, missing, 1), Prior(), 100_000)

pN = sample(sample_priorN(missing, missing, 10), Prior(), 100_000)

density(logistic.(pT[:intercept] .+ pT[:pc1] .+ pT[:pc2] .+ pT[:pc3] .+ pT[:pc4] .+ pT[:pc5] .+ pT[:pc6]))

density(logistic.(pN[:intercept] .+ pN[:pc1] .+ pN[:pc2] .+ pN[:pc3] .+ pN[:pc4] .+ pN[:pc5] .+ pN[:pc6]))