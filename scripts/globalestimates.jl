module CountOutputs

export nothing

using CSV, DataFrames, Chain

import StatsBase: denserank, mean

const ROOT = dirname(Base.active_project())

include("$ROOT/src/global.jl")
using .GlobalVariables

known = pop_known[:, [:atoll, :region, :species, :nbirds]]
glob = CSV.read("$ROOT/data/atoll_seabird_global_popestimates.csv", DataFrame)

DataFrames.transform!(pop_unknown, [:atoll, :species, :region] .=> denserank => x -> string("num_", x))

# Load data sets
p75, p80, p85 = CSV.read.("$ROOT/results/data/countpreds_0." .* ["75", "8", "85"] .* ".csv", DataFrame)
select!.([p75, p80, p85], Ref(Not(:upper, :lower)))

full75, full80, full85 = [vcat(known, df) for df in [p75, p80, p85]]

tot75, tot80, tot85 = @chain [full75, full80, full85] begin
    groupby.(_, :species)
    combine.(_, :nbirds => sum => identity)
    leftjoin.(_, Ref(glob), on=:species)
end

function calcratio(n, blmin, blmax, hbw, otero)
    ismissing(blmin) && ismissing(hbw) ? n / otero :
    ismissing(blmin) ? n / hbw :
    ismissing(blmax) ? n / blmin :
    n / mean([blmin, blmax])
end

tot_species.ratio = map(calcratio, eachrow(tot_species))

transform!.([tot75, tot80, tot85], [:nbirds, :birdlife_min, :birdlife_max, :HBW, :Otero] => ByRow(calcratio) => :ratio)

histogram((tot75.ratio .- tot85.ratio) .* 100, title="Pop on atoll ratio diff .75 - .85", label=:none)
xticks!(-1:3, string.(-1:3) .* "%")

CSV.write("$ROOT/results/data/allpopulations.csv", full80)
CSV.write("$ROOT/results/data/ratios.csv", tot80)

end
