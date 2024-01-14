module PCs

export nothing

using CSV, DataFrames, Statistics, Turing, Chain, StatsPlots, Random
import MultivariateStats as MS
import StatsBase as SB

const ROOT = dirname(Base.active_project())

const PC_NAMES = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"]

envs = CSV.read("$ROOT/data/envs_jicp.csv", DataFrame, missingstring="NA")

DataFrames.transform!(envs, :human_population => ByRow(x -> ifelse(ismissing(x), round(median(skipmissing(envs.human_population)), digits=0), x)) => identity)

DataFrames.transform!(envs, [:number_islets, :land_area_sqkm, :distance_nearest_atoll_km, :distance_nearest_high_island_km] .=> ByRow(x -> log(x)) => identity)
DataFrames.transform!(envs, [:lagoon_area_sqkm, :tropical_storms_50km, :hurricanes_50km, :human_population] .=> ByRow(x -> log1p(x)) => identity)

seabirds = @chain CSV.read("$ROOT/data/atoll_seabird_populations.csv", DataFrame) begin
    stack(_, Not(:atoll), variable_name=:species, value_name=:presence)
    DataFrames.transform(_, :presence => ByRow(x -> ismissing(x) ? 0.0 : 1.0) => identity)
end

df_features = envs[!, [collect(7:12)..., 15, collect(17:25)...]] # was 8:13, 
X_features = Matrix{Float64}(df_features)'

Random.seed!(42)

Z = MS.fit(SB.ZScoreTransform, X_features)
SB.transform!(Z, X_features)
# M12 = MS.fit(MS.PCA, X_features; maxoutdim=12) # We will work with the first 6 PCs
M = MS.fit(MS.PCA, X_features; maxoutdim=6)

envscores = @chain predict(M, X_features)' begin
    DataFrame(_, PC_NAMES)
    insertcols(_, :atoll => envs.atoll)
    insertcols(_, :region => envs.region)
    outerjoin(_, seabirds, on=:atoll)
end

CSV.write("$ROOT/results/data/jl_envscores.csv", envscores)

# Make a heatmap from PCA projections
proj = MS.projection(M)
features = names(df_features)

df_proj = @chain begin
    DataFrame(proj, PC_NAMES)
    insertcols(_, 1, :features => features)
end

plotfeatures = [
    "Number islets",
    "Annual precipitation",
    "Land area",
    "Lagoon area",
    "Tropical storms within 50km",
    "Hurricanes within 50km",
    "Distance nearest atoll",
    "Distance nearest high island",
    "Distance continent",
    "Human population",
    "Precipitation anomaly",
    "Total primary production Phyto",
    "Total chlorophyll",
    "Total phytoplankton",
    "Sea surface tempterature",
    "Wind speed",
]

cg = palette([:blue, :white, :white, :red], 25)

heatmap(PC_NAMES, collect(1:16), proj, size=(800, 600), color=cg, clim=(-0.75, 0.75))
yticks!(collect(1:16), plotfeatures)

foreach(ext -> savefig("$ROOT/results/$ext/pca.$ext"), ["svg", "png"])

end
