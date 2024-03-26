# This script is part of the project associated with
# Article: Atolls are globally significant sites for tropical seabirds
# Authors: Steibl S, Steiger S, Wegmann AS, Holmes ND, Young, HS, Carr P, Russell JC 
# Last edited: 2024-03-24

# Depends on count.jl, so make sure it is loaded
isdefined(Main, :CountModel) || include(joinpath(Main.ROOT, "julia", "scripts", "count.jl"))

"""
This module calculates summarises the results of the count pipeline by
- calculating global estimates
- running sensitivity analyses
- calculating nutrient input
"""
module Summaries

using Chain, DataFrames, CSV, StatsPlots, OrderedCollections, StatsBase

# Load custom modules
include(joinpath(Main.ROOT, "julia", "src", "utilities.jl"))
using .CustomUtilityFuns
include(joinpath(Main.ROOT, "julia", "src", "globalvars.jl"))
using .GlobalVariables

# --- GLOBAL ESTIMATES --- #

known = @chain pop_known begin
    transform(_, :nbirds => ByRow(x -> (median=x, lower=x, upper=x)) => AsTable)
    select(_, [:atoll, :region, :species, :median, :lower, :upper])
end

glob = CSV.read(joinpath(Main.ROOT, "data", "birdlife_hbw_globalestimates_$(Main.SUFFIX).csv"), DataFrame)

# Load data sets
p75, p80, p85 = CSV.read.(joinpath.(Main.ROOT, "results", "data", "countpreds_0." .* ["75", "8", "85"] .* "_default_$(Main.SUFFIX).csv"), DataFrame)
select!.([p75, p80, p85], Ref([:atoll, :region, :species, :median, :lower, :upper]))

full75, full80, full85 = [vcat(known, df) for df in [p75, p80, p85]]

tot75, tot80, tot85 = popsum.([full75, full80, full85])
leftjoin!.([tot75, tot80, tot85], Ref(glob), on=:species)

for n in [:median, :lower, :upper], df in [tot75, tot80, tot85]
    transform!(df, [:species, n, :birdlife_min, :birdlife_max, :HBW, :Otero] => ByRow(calcratio) => string("ratio_", n))
end
select!.([tot75, tot80, tot85], Ref(Not(:birdlife_min, :birdlife_max, :HBW, :Otero)))

histogram((tot75.ratio_median .- tot85.ratio_median) .* 100, title="Pop on atoll ratio diff .75 - .85", label=:none)
xticks!(-1:3, string.(-1:3) .* "%")

CSV.write(joinpath(Main.ROOT, "results", "data", "pred_and_obs_atolls_$(Main.SUFFIX).csv"), full80)

# Join posterior samples onto full data frame
full_samples = @chain known begin
    transform(_, :median => ByRow(x -> fill(x, Int(Main.CountModel.nsamples * Main.CountModel.nchains))) => :raw)
    select(_, :atoll, :region, :species, :raw, :median, :lower, :upper)
    vcat(Main.CountModel.preds_target["0.8"], _)
end

summary_count =
    let indices = [:species, :atoll]
        out = map(indices, [0.8, 0.95]) do index, prob
            out = map(i -> summariseby(i, index, full_samples; prob=prob), unique(full_samples[:, index]))
            DataFrame(out)
        end
        Dict(Pair.(string.(indices), out))
    end

foreach(k -> CSV.write(joinpath(Main.ROOT, "results", "data", "summary_count_$(k)wise_$(Main.priorsetting)_$(Main.SUFFIX).csv"), summary_count[k]), keys(summary_count))

# --- SENSITIVITY ANALYSIS --- # 

# for 0.8 cutoffs

dict_sensitivity = @chain begin
    map(Main.COUNT_PRIORSETTINGS) do ps # always run all priorsettings here
        Pair(ps, CSV.read(joinpath(Main.ROOT, "results", "data", "countpreds_0.8_$(ps)_$(Main.SUFFIX).csv"), DataFrame))
    end
    Dict(_)
end

select!.(values(dict_sensitivity), Ref([:atoll, :region, :species, :median, :lower, :upper]))
[dict_sensitivity[k] = vcat(known, v) for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))]
[dict_sensitivity[k] = popsum(v) for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))]
leftjoin!.(values(dict_sensitivity), Ref(glob), on=:species)

for n in [:median, :lower, :upper], df in values(dict_sensitivity)
    transform!(df, [:species, n, :birdlife_min, :birdlife_max, :HBW, :Otero] => ByRow(calcratio) => string("ratio_", n))
end

for (k, v) in Pair.(keys(dict_sensitivity), values(dict_sensitivity))
    dict_sensitivity[k] = insertcols(v, 1, :prior => k)
end

# Compare global predictions for different prior settings
df_sens = reduce(vcat, values(dict_sensitivity))

scatter(df_sens.ratio_median, df_sens.species, group=df_sens.prior, alpha=0.75, markershape=:x)
xlabel!("Global population on atolls")
xticks!(0:0.5:2, string.(Int64.(collect(0:0.5:2) .* 100), "%"))
yticks!(eachindex(unique(df_sens.species)), unique(df_sens.species), size=(600, 600), tickfontsize=7)

foreach(ext -> savefig(joinpath(Main.ROOT, "results", ext, "count", "sensitivity_count.$ext")), ["svg", "png"])

# --- NUTRIENT INPUT --- #

comb = leftjoin(atollinfo, select(full_samples, :atoll, :species, :raw), on=:atoll)
subset!(comb, :species => ByRow(x -> !ismissing(x))) # atolls with neither seabird data nor any predicted species
leftjoin!(comb, specinfo, on=:species)

const F_nc = 0.036 # average N content of seabird diet in g N g^-1
const F_pc = 0.006 # average P content of seabird diet in g P g^-1
const F_nv = 0.6   # proportion of volatilised excreted N
const F_ec = 6.5   # average energy content of seabird diet in kJ g^-1
const A_eff = 0.8  # kJ obtained by the bird per kJ consumed
const massratio = 17 / 14 # mass ratio of NH3 to N
const H2O_ratio = 0.7 # percent water wet weight
const C_ratio = 0.5 # carbon dry weight

# Adults
# Basic metabolic rate (kJ bird^-1 day^-1)
BMR = 2.3 .* comb.bodymass .^ 0.774

# Field metabolic rate (kJ bird^-1 day^-1)
AMR = 4 .* BMR

# Total daily N and P excreted (g N or g P bird^-1 day^-1)
adult_dN, adult_dP = [(AMR .* Fx) ./ (F_ec .* A_eff) for Fx in [F_nc, F_pc]]
adult_sN, adult_sP = [dx .* comb.days_at_colony .* comb.time_at_colony for dx in [adult_dN, adult_dP]]

# Extrapolate for colony size
adult_tN, adult_tP = [sx .* comb.raw for sx in [adult_sN, adult_sP]]

# Chicks
# Metabolic rate for chicks (kJ bird^-1 year^-1)
E_rear = 28.43 .* comb.fledgling_mass .^ 1.06

# Total daily N and P excreted (g N or g P bird^-1 day^-1)
chick_dN, chick_dP = [((E_rear .* Fx) ./ (F_ec .* A_eff)) .* (comb.prod_per_pair ./ 2) for Fx in [F_nc, F_pc]]

# Extrapolate for colony size
chick_tN, chick_tP = [sx .* comb.raw for sx in [chick_dN, chick_dP]]

# Add up adults and chicks and convert to kg
excretedN = (adult_tN .+ chick_tN) ./ 1000
excretedP = (adult_tP .+ chick_tP) ./ 1000

# NH3 emissions
F_hab = (comb.nestingtype .!= "burrow") .* 0.2 # re-adsorption of NH3 by nesting type
NH3emit = excretedN .* F_nv .* massratio .* F_hab

# Total bird biomass
biomass = comb.bodymass .* comb.raw

bird_C = biomass .* H2O_ratio .* C_ratio

# Summarise atoll-wise

df_raw_nutrients = DataFrame((; biomass, excretedN, excretedP, NH3emit, bird_C))
df_comb_nutr = [comb df_raw_nutrients]

summary_nutrient = let indices = [:atoll, :species]
    out = map(indices) do index
        vdf_nutrients_by_index = map(names(df_raw_nutrients)) do nutrient
            out = map(i -> summariseby(i, index, df_comb_nutr; target=nutrient), unique(df_comb_nutr[:, index]))
            select(DataFrame(out), index, Not(index) .=> identity => (x -> string(x, "_$nutrient")))
        end
        reduce(hcat, [select(vdf_nutrients_by_index[1], index), select.(vdf_nutrients_by_index[2:end], Ref(Not(index)))...])
    end
    Dict(Pair.(string.(indices), out))
end

foreach(k -> CSV.write(joinpath(Main.ROOT, "results", "data", "summary_nutrient_$(k)wise_$(Main.priorsetting)_$(Main.SUFFIX).csv"), summary_nutrient[k]), keys(summary_nutrient))

# Create global summaries for text
summary_global = map(nutrient -> summariseby(nothing, df_comb_nutr; target=nutrient), names(df_raw_nutrients))
global_dict = Dict(Pair.(names(df_raw_nutrients), summary_global))
push!(global_dict, Pair("count", summariseby(nothing, full_samples)))

textsummary = stringify(global_dict, prefix="FORMAT [lower, median, upper]\n\n")

# Write text-based summary
open(joinpath(Main.ROOT, "results", "data", "summary_textbased.txt"), "w") do file
    write(file, textsummary)
end

end