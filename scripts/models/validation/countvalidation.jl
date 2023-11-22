using CSV, Serialization, DataFrames

# Pathing
const ROOT = dirname(Base.active_project())

include("$ROOT/scripts/visualization/diagnosticplots.jl")
using .DiagnosticPlots

chain_narrow = deserialize("$ROOT/scripts/models/chains/chains_count_narrow.jls")
chain_wide = deserialize("$ROOT/scripts/models/chains/chains_count_wide.jls")
chain_default = deserialize("$ROOT/scripts/models/chains/chains_count_default.jls")

preds_narrow = CSV.read("$ROOT/data/countpreds_narrow.csv", DataFrame)
preds_wide = CSV.read("$ROOT/data/countpreds_wide.csv", DataFrame)
preds_default = CSV.read("$ROOT/data/countpreds_default.csv", DataFrame)

using StatsPlots, Chain

min_set_atoll = Set(preds_narrow.atoll)
min_set_species = Set(preds_narrow.species)

map([:atoll, :species], [min_set_atoll, min_set_species]) do var, set
    subset!.([preds_default, preds_narrow, preds_wide], var => ByRow(x -> x ∈ set))
end

long_default, long_narrow, long_wide = 
    @chain begin
    [preds_default, preds_narrow, preds_wide]
    unstack.(_, :species, :nbirds)
    transform.(_, All() .=> ByRow(x -> ismissing(x) ? 0 : x) => identity)
    select.(_, Ref(Not(:atoll, :region)))
    Matrix.(_)
end

diffs = long_default .- long_wide
heatmap(diffs)
extrema(diffs)
