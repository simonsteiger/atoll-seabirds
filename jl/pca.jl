using CSV, DataFrames, Statistics
import MultivariateStats as MS
import StatsBase as SB
using Match
# using MLDataUtils

envs = @chain CSV.read("data/envs_jicp.csv", DataFrame, missingstring="NA") begin
    DataFrames.transform(_, :human_population => ByRow(x -> ifelse(ismissing(x), round(median(skipmissing(envs.human_population)), digits=0), x)) => identity)
end

seabirds = @chain CSV.read("data/atoll_seabird_populations_11Mar.csv", DataFrame) begin
    stack(_, Not(:atoll), variable_name=:species, value_name=:presence)
    DataFrames.transform(_, :presence => ByRow(x -> ismissing(x) ? 0.0 : 1.0) => identity)
end

cond = CSV.read("data/seabird_filterconditions_03Jul.csv", DataFrame)

X_envs = Matrix{Float64}(envs[!, [collect(8:13)..., 16, collect(18:26)...]])'

Z = MS.fit(SB.ZScoreTransform, X_envs)
SB.transform!(Z, X_envs)
# M12 = MS.fit(MS.PCA, X_envs; maxoutdim=12) # We will work with the first 6 PCs
M = MS.fit(MS.PCA, X_envs; maxoutdim=6)

v_cond = Vector(fill("", nrow(cond)))

for i in eachindex(cond.filtercondition)
    v_cond[i] = @match cond.filtercondition[i] begin
        "none" => "\\w+"
        "EXCLUDE" => "\\d"
        _ => cond.filtercondition[i]
    end
end

cond.filtercondition = replace.(v_cond, r",\s" => "|")

envscores = @chain predict(M, X_envs)' begin
    DataFrame(_, [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6])
    insertcols(_, :atoll => envs.atoll)
    insertcols(_, :region => envs.region)
    outerjoin(_, seabirds, on=:atoll)
    outerjoin(_, cond, on=:species, matchmissing=:equal)
end

CSV.write("data/jl_envscores.csv", envscores)
