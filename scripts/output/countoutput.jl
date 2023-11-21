module CountOutputs

export nothing

using CSV, DataFrames, Chain

import StatsBase: denserank, mean

const ROOT = dirname(Base.active_project())

include("$ROOT/scripts/preprocessing/preprocess.jl")

known = Preprocess.pop_known[:, [:atoll, :region, :species, :nbirds]]
glob = CSV.read("$ROOT/data/atoll_seabird_global_popestimates.csv", DataFrame)

DataFrames.transform!(Preprocess.pop_unknown, [:atoll, :species, :region] .=> denserank => x -> string("num_", x))

# Load data sets
preds = @chain begin
    CSV.read("$ROOT/data/countpreds.csv", DataFrame)
    DataFrames.transform(_, :nbirds => ByRow(exp) => identity)
    rename(_, :atoll => :num_atoll, :species => :num_species, :region => :num_region)
end

helpjoin!(var) = leftjoin!(preds, unique(Preprocess.pop_unknown[:, [var, "num_$var"]], var), on="num_$var")

map(helpjoin!, ["atoll", "species", "region"])

select!(preds, [:atoll, :species, :region, :nbirds], r"^(?!num_).+")

full = vcat(known, preds)

tot_species = @chain full begin
    groupby(_, :species)
    combine(_, :nbirds => sum => identity)
    leftjoin(_, glob, on=:species)
end

function calcratio(row)
    ismissing(row.birdlife_min) && ismissing(row.HBW) ? row.nbirds / row.Callaghan :
    ismissing(row.birdlife_min) ? row.nbirds / row.HBW :
    ismissing(row.birdlife_max) ? row.nbirds / row.birdlife_min :
    row.nbirds / mean([row.birdlife_min, row.birdlife_max])
end

tot_species.ratio = map(calcratio, eachrow(tot_species))

end
