# Include local modules
include("/Users/simonsteiger/Desktop/other/atoll-seabirds/jl/aux_prep.jl")
include("/Users/simonsteiger/Desktop/other/atoll-seabirds/jl/aux_model.jl")

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

Random.seed!(0);

# Set prior degrees of freedom for Student T distribution
const PRIOR_DF = 3

struct result
    chain
    threshold
    diagnostics
end

# Import the data
envscores = CSV.read("data/jl_envscores.csv", DataFrame)

envs_known = subset(envscores, :presence => ByRow(x -> !ismissing(x)))
envs_unknown = subset(envscores, :presence => ByRow(x -> ismissing(x)))

# Delete all species with known population from data
envs_known = subset(envs_known, [:filtercondition, :region] => ByRow((x, y) -> occursin(Regex(x), y)))

# Discard unused columns
select!(envs_known, Not([:region, :filtercondition]))

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

# all_results_N = Dict{String,result}()
all_results = Dict{String,result}()

[all_results[k] = wrap_model(k, logistic_regression, PRIOR_DF) for k in collect(keys(trainset))]

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

# TODO
# Check if there is a correlation between low majority class N and prediction performance
# Can we fit a single model with indices for each species?

