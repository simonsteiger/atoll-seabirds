using Turing, CSV, DataFrames, Chain

const ROOT = dirname(Base.active_project())

preds = CSV.read("$ROOT/data/allpopulations.csv", DataFrame)
envs = CSV.read("$ROOT/data/seabird_atolls_envs.csv", DataFrame)

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

preds.Δ = rand(Uniform(15, 150), nrow(preds))

preds.nbirds = Int64.(round.(preds.nbirds))

distance = @chain preds [fill(Weibull(8, r.Δ), r.nbirds) for r in eachrow(_)] [rand.(W) for W in _] reduce(vcat, _)

direction = @chain preds.nbirds [fill(Uniform(0, 360), n) for n in _] [rand.(U) for U in _] reduce(vcat, _)

atolls = reduce(vcat, [fill(preds.atoll[i], preds.nbirds[i]) for i in eachindex(preds.nbirds)])

out = DataFrame([atolls, distance, direction], [:atoll, :distance, :direction])

DataFrames.transform!(out, [:direction, :distance] => ByRow((dir, dis) -> (addlatkm=cos((dir * π) / 180) * dis, addlongkm=sin((dir * π) / 180) * dis)) => AsTable)

leftjoin!(out, select(envs, :atoll, :lat, :long), on=:atoll)

out.long = [long < 0 ? long + 360 : long for long in out.long]

out.addlatdeg = vec(reduce(vcat, km2d.(out.addlatkm, out.lat)))
out.addlongdeg = vec(reduce(vcat, km2d.(out.addlongkm, out.long)))

select!(out, :lat, :long, :addlatdeg => :addlat, :addlongdeg => :addlong)

transform!(out, [[:lat, :addlat], [:long, :addlong]] .=> ByRow((x,y) -> round(x + y, digits=1)) .=> [:adjlat, :adjlong])
select!(out, r"adj")

subset_size = round(Int, 0.1 * nrow(out))
random_rows = sample(eachindex(out.adjlat), subset_size, replace = false)

sub_out = out[random_rows, :]

#out = @chain out begin
#    groupby(_, [:adjlat, :adjlong])
#    combine(_, nrow)
#    transform(_, :nrow => ByRow(log) => :dens)
#    select(_, Not(:nrow))
#end

using RCall

@rput out

R"""
library(fst)
library(here)

write_fst(out, here("results/data/foraging.fst"))
"""
