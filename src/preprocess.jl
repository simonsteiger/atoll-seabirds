module Preprocess

export train, test, train_label, test_label, groupatoll, envs_known,
       all_atoll, all_species, all_presence, all_pc, all_nesting, 
       distM, distM_known, distM_unknown

# Working with tabular data
using CSV, DataFrames, Chain

using StatsBase

using MLDataUtils: shuffleobs, stratifiedobs, oversample, rescale!
using Resample

using Random, StableRNGs

include("upsample.jl")
include("tune.jl")

Random.seed!(0);

const PRIOR_TEST = 0.5
const PRIOR_DF = 3

# Import the data
envscores = CSV.read("data/jl_envscores.csv", DataFrame)

# Import the data
distances = CSV.read("data/pairwise_atoll_distance.csv", DataFrame)
distM = Matrix{Float64}(distances[:, 2:end])
distM ./= 1000

# Split distM into distM_known and distM_unknoqn
ispresencemissing = @chain envscores begin
    unique(_, :atoll)
    getproperty(_, :presence)
    map(x -> ismissing(x), _)
end

distM_known = distM[.!ispresencemissing, .!ispresencemissing]
distM_unknown = distM[ispresencemissing, ispresencemissing]

# Add data about nesting type
specinfo = @chain begin
    CSV.read("data/seabird_filterconditions_13Sep.csv", DataFrame)
    select(_, [:species, :nestingtype])
end

# Convert "presence" to numeric values
# envscores[!, :presence] = [r.presence == true ? 1.0 : 0.0 for r in eachrow(envscores)]

envs_known = subset(envscores, :presence => ByRow(x -> !ismissing(x)))
envs_unknown = subset(envscores, :presence => ByRow(x -> ismissing(x)))

# Delete all species with known population from data
# RETHINK Do we want to keep known populations to inform other groups?
# envs_known = subset(envs_known, [:filtercondition, :region] => ByRow((x, y) -> occursin(Regex(x), y)))

# Add nestingtype to envscores
leftjoin!(envs_known, specinfo, on=:species)

# Discard unused columns
envs_model = select(envs_known, Not([:region, :filtercondition, :nestingtype]))

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
# trainset = Dict{String,DataFrame}()
# testset = Dict{String,DataFrame}()
# 
# [(trainset[s], testset[s]) = split_data(envs_model, target, s; at=0.5) for s in unique(envs_model.species)]
# 
# for f in numerics, k in keys(trainset)
#     μ, σ = rescale!(trainset[k][!, f]; obsdim=1)
#     rescale!(testset[k][!, f], μ, σ; obsdim=1)
# end
# 
# trainset_up = Dict{String,DataFrame}()
# 
# rng = StableRNG(1)
# 
# [trainset_up[k] = our_smote(rng, trainset[k]) for k in keys(trainset)];
# 
# # Dicts for train, test, train_label, test_label
# train = Dict{String,Matrix}()
# test = Dict{String,Matrix}()
# train_label = Dict{String,Vector}()
# test_label = Dict{String,Vector}()
# groupatoll = Dict{String, Vector}()
# 
# for k in keys(trainset)
#     train[k] = Matrix(trainset[k][:, features])
#     test[k] = Matrix(testset[k][:, features])
#     train_label[k] = trainset[k][:, target]
#     test_label[k] = testset[k][:, target]
#     groupatoll[k] = denserank(trainset[k][:, :atoll])
# end

all_atoll = Int64.(denserank(envs_known.atoll))
all_species = Int64.(denserank(envs_known.species))
all_presence = Float64.(envs_known.presence)
all_pc = Matrix(envs_known[:, features])
all_nesting = Int64.(denserank(envs_known.nestingtype))

end