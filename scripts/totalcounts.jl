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
p75, p80, p85 = CSV.read.("$ROOT/results/data/countpreds_0." .* ["75", "8", "85"] .* "_default.csv", DataFrame)
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

CSV.write("$ROOT/results/data/allpopulations.csv", full80)
CSV.write("$ROOT/results/data/ratios.csv", tot80)

# --- SENSITIVITY ANALYSIS --- # 

# for 0.8 cutoffs

dict_sensitivity = @chain begin
    map(["default", "mean", "global", "narrow", "wide"]) do prior
        Pair(prior, CSV.read("$ROOT/results/data/countpreds_0.8_$prior.csv", DataFrame))
    end
    Dict(_)
end

select!.(values(dict_sensitivity), Ref([:atoll, :region, :species, :median, :lower, :upper]))
rename!.(values(dict_sensitivity), :median => :nbirds)
[dict_sensitivity[k] = vcat(known, v) for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))]
[dict_sensitivity[k] = popsum(v) for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))]
leftjoin!.(values(dict_sensitivity), Ref(glob), on=:species)

for n in [:nbirds, :lower, :upper], df in values(dict_sensitivity)
    transform!(df, [n, :birdlife_min, :birdlife_max, :HBW, :Otero] => ByRow(calcratio) => string("ratio_", n))
end

for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))
    dict_sensitivity[k] = insertcols(v, 1, :prior => k)
end

df_sens = reduce(vcat, values(dict_sensitivity))

scatter(df_sens.ratio_nbirds, df_sens.species, group=df_sens.prior, alpha=0.75, markershape=:x)
xlabel!("Global population on atolls")
xticks!(0:0.5:2, string.(Int64.(collect(0:0.5:2) .* 100), "%"))
yticks!(eachindex(unique(df_sens.species)), unique(df_sens.species), size=(600,600), tickfontsize=7)
savefig("$ROOT/results/svg/sensitivity_count.svg")

end
