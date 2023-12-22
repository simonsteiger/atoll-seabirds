module GlobalVariables

export envs_known,
       envs_unknown,
       pop_known,
       pop_unknown

# Working with tabular data
using CSV, DataFrames, Chain

import StatsBase: denserank

ARGS2 = isempty(ARGS) ? "default" : ARGS[3]

ROOT = dirname(Base.active_project())

# Import the data
envscores = CSV.read("$ROOT/results/data/jl_envscores.csv", DataFrame)
pop = CSV.read("$ROOT/data/atoll_seabird_populations.csv", DataFrame)

# Add data about nesting type
specinfo = @chain begin
    CSV.read("$ROOT/data/seabird_filterconditions.csv", DataFrame)
    select(_, [:species, :nestingtype])
end

envs_known = subset(envscores, :presence => ByRow(x -> !ismissing(x)))
envs_unknown = subset(envscores, :presence => ByRow(x -> ismissing(x)))

# Dict for back-and-forth translating between index and string notation
odict_atoll_unknown = sort(Dict(Pair.(denserank(envs_unknown.atoll), envs_unknown.atoll)))

# Add nestingtype to envscores
leftjoin!(envs_known, specinfo, on=:species)

pop_known = @chain pop begin
    DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    stack(_, Not(:atoll), variable_name=:species, value_name=:nbirds)
    subset(_, :nbirds => ByRow(x -> x != -1))
    subset(_, :nbirds => ByRow(x -> x != 0))
    leftjoin(_, envs_known, on=[:atoll, :species])
    select(_, Not(:presence))
end

# Dict for back-and-forth translating between index and string notation
odict_species_known = sort(Dict(Pair.(denserank(pop_known.species), pop_known.species)))

# Create pop_unknown
pop_unknown = @chain pop begin
    DataFrames.transform(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    stack(_, Not(:atoll), variable_name=:species, value_name=:nbirds)
    subset(_, :nbirds => ByRow(x -> x == -1))
    leftjoin(_, dropmissing!(envscores, :species), on=[:atoll, :species])
    select(_, Not(:presence))
end



# Check if prediction data frame for a given prior setting exists
preds_exist = isfile("$ROOT/results/data/presencepreds.csv")

# Create pred data frame if this is the case
if preds_exist
    preds = @chain begin
        CSV.read("$ROOT/results/data/presencepreds.csv", DataFrame)
        transform(_, :atoll => ByRow(x -> haskey(odict_atoll_unknown, x) ? odict_atoll_unknown[x] : x) => identity)
        transform(_, :species => ByRow(x -> haskey(odict_species_known, x) ? odict_species_known[x] : x) => identity)
        select(_, :atoll, :mean => :nbirds, r"PC", :species) # rename to nbirds for vcat below, is actually prob present (ppres)
        leftjoin(_, select(envs_unknown, r"PC|region|atoll"), on=:atoll)
    end
else
    preds = nothing
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