module Preprocess

export envs_known,
    envs_unknown,
    pop_known,
    pop_unknown

# Working with tabular data
using CSV, DataFrames, Chain

import StatsBase: denserank

ARGS2 = isempty(ARGS) ? "default" : ARGS[2]

ROOT = dirname(Base.active_project())

# Import the data
envscores = CSV.read("$ROOT/data/jl_envscores.csv", DataFrame)

# Split distM into distM_known and distM_unknoqn
ispresencemissing = @chain envscores begin
    unique(_, :atoll)
    getproperty(_, :presence)
    map(x -> ismissing(x), _)
end

# Add data about nesting type
specinfo = @chain begin
    CSV.read("$ROOT/data/seabird_filterconditions_13Sep.csv", DataFrame)
    select(_, [:species, :nestingtype])
end

# Convert "presence" to numeric values
# envscores[!, :presence] = [r.presence == true ? 1.0 : 0.0 for r in eachrow(envscores)]

envs_known = subset(envscores, :presence => ByRow(x -> !ismissing(x)))
envs_unknown = subset(envscores, :presence => ByRow(x -> ismissing(x)))

# Add nestingtype to envscores
leftjoin!(envs_known, specinfo, on=:species)

pop_known = @chain begin
    CSV.read("$ROOT/data/atoll_seabird_populations_29Jul.csv", DataFrame)
    DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    subset(_, All() .=> ByRow(x -> x != -1))
    stack(_, Not(:atoll), variable_name=:species, value_name=:nbirds)
    subset(_, :nbirds => ByRow(x -> x != 0))
    leftjoin(_, envs_known, on=[:atoll, :species])
    select(_, Not(:presence))
end

# Check if prediction data frame for a given prior setting exists
preds_exist = isfile("$ROOT/data/presencepreds_$ARGS2.csv")

# Create pred data frame if this is the case
if preds_exist
    @info "Using $ARGS2 predictions."
    preds = @chain begin
        CSV.read("$ROOT/data/presencepreds_$ARGS2.csv", DataFrame)
        stack(_)
        select(_, :atoll, :variable => :species, :value => :nbirds)
        leftjoin(_, select(envs_unknown, r"PC|region|atoll"), on=:atoll)
    end
else
    preds = nothing
end

# Create pop_unknown
pop_unknown = @chain begin
    CSV.read("$ROOT/data/atoll_seabird_populations_29Jul.csv", DataFrame)
    DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    stack(_, Not(:atoll), variable_name=:species, value_name=:nbirds)
    subset(_, :nbirds => ByRow(x -> x == -1))
    leftjoin(_, dropmissing!(envscores, :species), on=[:atoll, :species])
    select(_, Not(:presence))
end

# Add predictions to pop_unknown if they exist for the given prior setting
if preds_exist
    # DataFrame to join real species name back onto numeric values
    df_species_nestingtype = unique(envs_known[:, [:species, :nestingtype]], :species)
    # Add predictions to pred_unknown and add full species names back
    pop_unknown = @chain pop_unknown begin
        vcat(_, preds)
        subset(_, :species => ByRow(x -> x ∈ pop_known.species))
        select(_, :atoll, :region, :species, :nbirds => :ppres, Cols(contains("PC")))
        transform(_, :ppres => ByRow(x -> x == -1 ? 1.0 : x) => identity)
        leftjoin(_, df_species_nestingtype, on=:species)
    end
end

end