module Preprocess

export envs_known,
       envs_unknown,
       pop_known,
       pop_unknown

# Working with tabular data
using CSV, DataFrames, Chain

import StatsBase: denserank

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
    select(_, Not(:presence))
end

preds = @chain begin
    CSV.read("../../data/newpreds.csv", DataFrame)
    stack(_)
    select(_, :atoll, :variable => :species, :value => :nbirds)
    leftjoin(_, select(envs_unknown, r"PC|region|atoll"), on=:atoll)
end

df_species_nestingtype = unique(envs_known[:, [:species, :nestingtype]], :species)

pop_unknown = @chain begin
    CSV.read("../../data/atoll_seabird_populations_29Jul.csv", DataFrame)
    DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    stack(_, Not(:atoll), variable_name=:species, value_name=:nbirds)
    subset(_, :nbirds => ByRow(x -> x == -1))
    leftjoin(_, dropmissing!(envscores, :species), on=[:atoll, :species])
    select(_, Not(:presence))
    vcat(_, preds)
    subset(_, :species => ByRow(x -> x ∈ pop_known.species))
    select(_, :atoll, :region, :species, :nbirds => :ppres, Cols(contains("PC")))
    transform(_, :ppres => ByRow(x -> x == -1 ? 1.0 : x) => identity)
    leftjoin(_, df_species_nestingtype, on=:species)
end

end