using CSV, DataFrames, Statistics, Turing, Chain, StatsPlots
import MultivariateStats as MS
import StatsBase as SB
using Match

include("../../src/standardise.jl")

const PC_NAMES = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"]

envs = CSV.read("../../data/envs_jicp.csv", DataFrame, missingstring="NA")

DataFrames.transform!(envs, :human_population => ByRow(x -> ifelse(ismissing(x), round(median(skipmissing(envs.human_population)), digits=0), x)) => identity)

seabirds = @chain CSV.read("../../data/atoll_seabird_populations_11Mar.csv", DataFrame) begin
    stack(_, Not(:atoll), variable_name=:species, value_name=:presence)
    DataFrames.transform(_, :presence => ByRow(x -> ismissing(x) ? 0.0 : 1.0) => identity)
end

cond = CSV.read("../../data/seabird_filterconditions_03Jul.csv", DataFrame)

df_features = envs[!, [collect(8:13)..., 16, collect(18:26)...]]
X_features = Matrix{Float64}(df_features)'

Z = MS.fit(SB.ZScoreTransform, X_features)
SB.transform!(Z, X_features)
# M12 = MS.fit(MS.PCA, X_features; maxoutdim=12) # We will work with the first 6 PCs
M = MS.fit(MS.PCA, X_features; maxoutdim=6)

v_cond = Vector(fill("", nrow(cond)))

for i in eachindex(cond.filtercondition)
    v_cond[i] = @match cond.filtercondition[i] begin
        "none"      => "\\w+"
        "EXCLUDE"   => "\\d"
        _           => cond.filtercondition[i]
    end
end

cond.filtercondition = replace.(v_cond, r",\s" => "|")

envscores = @chain predict(M, X_features)' begin
    [standardise(col) for col in eachslice(_, dims=2)]
    DataFrame(_, PC_NAMES)
    insertcols(_, :atoll => envs.atoll)
    insertcols(_, :region => envs.region)
    outerjoin(_, seabirds, on=:atoll)
    outerjoin(_, cond, on=:species, matchmissing=:equal)
end

# CSV.write("data/jl_envscores.csv", envscores)

# Make a heatmap from PCA projections
proj = MS.projection(M)
features = names(df_features)

df_proj = @chain begin
    DataFrame(proj, PC_NAMES)
    insertcols(_, 1, :features => features)
end

c = palette(:viridis, 50)

cg = cgrad([c[1], :grey80, :grey80, c[50]])

heatmap(PC_NAMES, collect(1:16), proj, size=(800, 600), color=cg)
yticks!(collect(1:16), features)
