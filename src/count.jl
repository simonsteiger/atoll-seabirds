module CountVariables

export atoll, region, species, nesting, species_in_nesting, PC, nbirds, nspecies, ppres, oos_lims
       
include("global.jl")
using .GlobalVariables

using DataFrames, Chain, StatsBase, CSV, OrderedCollections

const FEATURES = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]

ROOT = dirname(Base.active_project())

# Import the outofsample validation data
df_oos = @chain "$ROOT/data/atoll_seabird_populations_outofsample-validation.csv" begin
    # Three columns formatted as date in original file, make sure to read as String
    CSV.read(_, DataFrame, types=Union{String, Missing})
    stack(_, Not(:atoll), variable_name=:species)
    dropmissing(_, :value)
    transform(_, :value => ByRow(x -> split(string(x), "-")) => :substrings)
    # CSV.read saves memory by encoding some strings as String7, but push!(...) requires String15 to work correctly (not sure why)
    transform(_, :substrings => ByRow(x -> any(occursin.(">", x)) ? push!(String15.(x), "Inf") : x) => identity)
    transform(_, :substrings => ByRow(x -> replace.(x, ">" => "")) => identity)
    transform(_, :substrings => ByRow(x -> parse.(Float64, x)) => :lims)
    leftjoin(_, unique(vcat(envs_known[:, [:atoll, FEATURES...]], envs_unknown[:, [:atoll, FEATURES...]])), on=:atoll)
    leftjoin(_, unique(pop_unknown[:, [:atoll, :region]]), on=:atoll)
    leftjoin(_, unique(pop_known[:, [:species, :nestingtype]], :species), on=:species)
end

oos_lims = df_oos.lims

# Repeat for population / count model

nspecies = @chain begin
    unique(pop_known, :species)
    groupby(_, :nestingtype)
    combine(_, nrow)
    getproperty(_, :nrow)
    (burrow=_[1], ground=_[2], vegetation=_[3])
end

# Atolls. ...
num_atoll_known, num_atoll_unknown = [Int64.(denserank(df.atoll)) for df in [pop_known, pop_unknown]]
atoll = (
    known = (num=num_atoll_known, str=pop_known.atoll),
    unknown = (num=num_atoll_unknown, str=pop_unknown.atoll),
    validation = (num=Int64.(denserank(df_oos.atoll)), str=df_oos.atoll)
)

# Region
odict_region = sort(Dict(Pair.(pop_known.region, Int64.(denserank(pop_known.region)))))

num_region_known, num_region_unknown, num_region_oos = 
    map([pop_known, pop_unknown, df_oos]) do df
        [haskey(odict_region, x) ? odict_region[x] : x for x in df.region]
    end

region = (
    known = (num=num_region_known, str=pop_known.region),
    unknown = (num=num_region_unknown, str=pop_unknown.region),
    validation = (num=Int64.(denserank(df_oos.region)), str=df_oos.region)
)

# Species
odict_species = sort(Dict(Pair.(pop_known.species, Int64.(denserank(pop_known.species)))))

num_species_known, num_species_unknown, num_species_oos = 
    map([pop_known, pop_unknown, df_oos]) do df
        [haskey(odict_species, x) ? odict_species[x] : x for x in df.species]
    end

species = (
    known = (num=num_species_known, str=pop_known.species),
    unknown = (num=num_species_unknown, str=pop_unknown.species),
    validation = (num=Int64.(denserank(df_oos.species)), str=df_oos.species)
)

nbirds = Float64.(pop_known.nbirds)
ppres = Float64.(pop_unknown.ppres)


PC_known = Matrix{Float64}(pop_known[:, FEATURES])
PC_unknown = Matrix{Float64}(pop_unknown[:, FEATURES])
PC_oos = Matrix{Float64}(df_oos[:, FEATURES])

PC = (
    known = PC_known,
    unknown = PC_unknown,
    validation = PC_oos,
)

# Nestingtype
odict_nesting = sort(Dict(Pair.(pop_known.nestingtype, Int64.(denserank(pop_known.nestingtype)))))

num_nesting_known, num_nesting_unknown, num_nesting_oos = 
    map([pop_known, pop_unknown, df_oos]) do df
        [haskey(odict_nesting, x) ? odict_nesting[x] : x for x in df.nestingtype]
    end

function unique_pop_var(df, unique_vars, pull)
    @chain df begin
        unique(_, unique_vars)
        DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
        getproperty(_, pull)
        denserank(_)
        map(Int64, _)
    end
end

unique_nesting_known = unique_pop_var(pop_known, [:nestingtype, :species], :nestingtype)
unique_species_in_nesting_known = unique_pop_var(pop_known, [:nestingtype, :species], :ne_sp)

nesting = (
    known = (num=num_nesting_known, str=pop_known.nestingtype),
    unknown = (num=num_nesting_unknown, str=pop_unknown.nestingtype),
    validation = (num=num_nesting_oos, str=df_oos.nestingtype),
    levels = unique_nesting_known,
)

# Species within nesting

dict_species_within_nesting = @chain pop_known begin
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    select(_, :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
    Dict(Pair.(pop_known.species, _))
end

num_species_in_nesting_known, num_species_in_nesting_unknown, num_species_in_nesting_oos =
    map([pop_known, pop_unknown, df_oos]) do df
        [haskey(dict_species_within_nesting, x) ? dict_species_within_nesting[x] : x for x in df.species]
    end

species_in_nesting = (
    known = (num=num_species_in_nesting_known),
    unknown = (num=num_species_in_nesting_unknown),
    validation = (num=num_species_in_nesting_oos),
    levels=unique_species_in_nesting_known
)

end