module GlobalEstimates

export nothing

using CSV, DataFrames, Chain
import StatsBase: denserank, mean
using StatsPlots

const ROOT = dirname(Base.active_project())

include("$ROOT/src/global.jl")
using .GlobalVariables

known = @chain pop_known begin
    select(_, [:atoll, :region, :species, :nbirds])
    transform(_, :nbirds => ByRow(x -> (lower=x, upper=x)) => AsTable)
end

glob = CSV.read("$ROOT/data/atoll_seabird_global_popestimates.csv", DataFrame)

DataFrames.transform!(pop_unknown, [:atoll, :species, :region] .=> denserank => x -> string("num_", x))

# Load data sets
p75, p80, p85 = CSV.read.("$ROOT/results/data/countpreds_0." .* ["75", "8", "85"] .* "_test.csv", DataFrame)
select!.([p75, p80, p85], Ref([:atoll, :region, :species, :median, :lower, :upper]))
rename!.([p75, p80, p85], :median => :nbirds)

full75, full80, full85 = [vcat(known, df) for df in [p75, p80, p85]]

function popsum(df)
    out = @chain df begin
        groupby(_, :species)
        combine(_, [:nbirds, :lower, :upper] .=> sum .=> identity)
    end
    return out
end

tot75, tot80, tot85 = popsum.([full75, full80, full85])
leftjoin!.([tot75, tot80, tot85], Ref(glob), on=:species)

function calcratio(n, blmin, blmax, hbw, otero)
    ismissing(blmin) && ismissing(hbw) ? n / otero :
    ismissing(blmin) ? n / hbw :
    ismissing(blmax) ? n / blmin :
    n / mean([blmin, blmax])
end

for n in [:nbirds, :lower, :upper], df in [tot75, tot80, tot85]
    transform!(df, [n, :birdlife_min, :birdlife_max, :HBW, :Otero] => ByRow(calcratio) => string("ratio_", n))
end
select!.([tot75, tot80, tot85], Ref(Not(:birdlife_min, :birdlife_max, :HBW, :Otero)))

histogram((tot75.ratio_nbirds .- tot85.ratio_nbirds) .* 100, title="Pop on atoll ratio diff .75 - .85", label=:none)
xticks!(-1:3, string.(-1:3) .* "%")

CSV.write("$ROOT/results/data/allpopulations_test.csv", full80)
CSV.write("$ROOT/results/data/ratios_test.csv", tot80)

end
