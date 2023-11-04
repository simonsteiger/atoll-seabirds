module Preprocess

export envs_known,
       envs_unknown,
       pop_known

# Working with tabular data
using CSV, DataFrames, Chain

using StatsBase

using MLDataUtils: shuffleobs, stratifiedobs, oversample, rescale!
using Resample

using Random, StableRNGs

include("../../src/upsample.jl")
include("../../src/tune.jl")

Random.seed!(0);

# Import the data
envscores = CSV.read("../../data/jl_envscores.csv", DataFrame)

# Split distM into distM_known and distM_unknoqn
ispresencemissing = @chain envscores begin
    unique(_, :atoll)
    getproperty(_, :presence)
    map(x -> ismissing(x), _)
end

# Add data about nesting type
specinfo = @chain begin
    CSV.read("../../data/seabird_filterconditions_13Sep.csv", DataFrame)
    select(_, [:species, :nestingtype])
end

# Convert "presence" to numeric values
# envscores[!, :presence] = [r.presence == true ? 1.0 : 0.0 for r in eachrow(envscores)]

envs_known = subset(envscores, :presence => ByRow(x -> !ismissing(x)))
envs_unknown = subset(envscores, :presence => ByRow(x -> ismissing(x)))

# Add nestingtype to envscores
leftjoin!(envs_known, specinfo, on=:species)

pop_known = @chain begin
    CSV.read("../../data/atoll_seabird_populations_29Jul.csv", DataFrame)
    DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    subset(_, All() .=> ByRow(x -> x != -1))
    stack(_, Not(:atoll), variable_name=:species, value_name=:nbirds)
    subset(_, :nbirds => ByRow(x -> x != 0))
    leftjoin(_, envs_known, on=[:atoll, :species])
    select(_, Not(:presence, :filtercondition))
end

end