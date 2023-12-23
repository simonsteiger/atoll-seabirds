using Turing, CSV, DataFrames, Chain

const ROOT = dirname(Base.active_project())

preds = CSV.read("$ROOT/data/allpopulations.csv", DataFrame)
envs = CSV.read("$ROOT/data/seabird_atolls_envs.csv", DataFrame)
raw_foraging = CSV.read("$ROOT/data/seabird_atolls_foraging.csv", DataFrame)

foraging = copy(raw_foraging)

select!(foraging, :species, r"maximum")

function addspecies(df, target, reference)
    return vcat(df, transform(df[df.species.==reference, :], :species => (_ -> target) => identity))
end

foraging = @chain foraging begin
    groupby(_, :species)
    combine(_, r"max" => (x -> (mean=mean(x), std=std(x))) => AsTable)
    transform(_, :species => ByRow(x -> occursin("Sternula_sp", x) ? "Sternula_saundersi" : x) => identity)
    transform(_, :species => ByRow(x -> occursin("Hydrobates", x) ? "Hydrobates_tristrami" : x) => identity)
    addspecies(_, "Nesofregetta_fuliginosa", "Hydrobates_tristrami")
    addspecies(_, "Onychoprion_lunatus", "Onychoprion_fuscatus")
    transform(_, :species => ByRow(x -> split(x, "_")[1]) => :genus)
    groupby(_, :genus)
    transform(_, [:mean, :std] => ((μ, σ) -> [isnan(s) ? std(μ) : s for s in σ]) => :std)
    addspecies(_, "Pterodroma_hypoleuca", "Pterodroma_alba")
    groupby(_, :genus)
    transform(_, [:species, :mean, :std] => ((sp, μ, σ) -> [occursin("hypoleuca", sp[i]) ? (mean=mean(μ), std=std(μ)) : (mean=μ[i], std=σ[i]) for i in eachindex(sp)]) => AsTable)
    transform(_, :genus => ByRow(x -> occursin(r"Hydroprogne|Thalasseus", x) ? "Hydroprogne/Thalasseus" : x) => identity)
    groupby(_, :genus)
    transform(_, [:species, :mean, :std] => ((sp, μ, σ) -> [occursin("caspia", sp[i]) ? (mean=μ[i], std=std(μ)) : (mean=μ[i], std=σ[i]) for i in eachindex(sp)]) => AsTable)
    select(_, Not(:genus))
end


function rdistearth(x1::AbstractArray{Float64}, x2::AbstractArray{Float64}, R=6378.388)
    coslat1 = @. cos((x1[:, 2] * π) / 180)
    sinlat1 = @. sin((x1[:, 2] * π) / 180)
    coslon1 = @. cos((x1[:, 1] * π) / 180)
    sinlon1 = @. sin((x1[:, 1] * π) / 180)

    coslat2 = @. cos((x2[:, 2] * π) / 180)
    sinlat2 = @. sin((x2[:, 2] * π) / 180)
    coslon2 = @. cos((x2[:, 1] * π) / 180)
    sinlon2 = @. sin((x2[:, 1] * π) / 180)

    pp = [coslat1 .* coslon1 coslat1 .* sinlon1 sinlat1] * [coslat2 .* coslon2 coslat2 .* sinlon2 sinlat2]'
    gt1 = pp .> 1
    pp[gt1] = 1 .* sign.(pp[gt1])

    return R .* acos.(pp)
end

km2d(km, baselat) = km ./ rdistearth([0.0 baselat], [1.0 baselat])

km2d.([-1.4, -1.2], [-2.3, -2.2])

preds.nbirds = Int64.(round.(preds.nbirds))
leftjoin!(preds, foraging, on=:species)

distance = @chain preds begin
    [fill(truncated(Normal(r.mean, r.std)), r.nbirds) for r in eachrow(_)]
    [rand.(W) for W in _]
    reduce(vcat, _)
end

direction = @chain preds.nbirds begin
    [fill(Uniform(0, 360), n) for n in _]
    [rand.(U) for U in _]
    reduce(vcat, _)
end

atolls = reduce(vcat, [fill(preds.atoll[i], preds.nbirds[i]) for i in eachindex(preds.nbirds)])

out = DataFrame([atolls, distance, direction], [:atoll, :distance, :direction])

DataFrames.transform!(out, [:direction, :distance] => ByRow((dir, dis) -> (addlatkm=cos((dir * π) / 180) * dis, addlongkm=sin((dir * π) / 180) * dis)) => AsTable)

leftjoin!(out, select(envs, :atoll, :lat, :long), on=:atoll)

out.long = [long < 0 ? long + 360 : long for long in out.long]

out.addlatdeg = vec(reduce(vcat, km2d.(out.addlatkm, out.lat)))
out.addlongdeg = vec(reduce(vcat, km2d.(out.addlongkm, out.long)))

# TODO subset to Cocos only and check if result is oval

select!(out, :lat, :long, :addlatdeg => :addlat, :addlongdeg => :addlong)

transform!(out, [[:lat, :addlat], [:long, :addlong]] .=> ByRow((x, y) -> round(x + y, digits=1)) .=> [:adjlat, :adjlong])
select!(out, r"adj")

using RCall

@rput out

R"""
library(fst)
library(here)

write_fst(out, here("results/data/foraging.fst"))
"""
