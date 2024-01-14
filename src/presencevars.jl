"This module creates objects necessary for fitting the presence model."
module PresenceVariables

export atoll, region, species, nesting, species_in_nesting, presence, PC, nspecies

include("globalvars.jl")
using .GlobalVariables

using DataFrames, Chain, StatsBase

const FEATURES = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]

nspecies = @chain begin
    unique(pop_known, :species)
    groupby(_, :nestingtype)
    combine(_, nrow)
    getproperty(_, :nrow)
    (burrow=_[1], ground=_[2], vegetation=_[3])
end

df_species_unknown = @chain begin
    unique(envs_known, :species)
    select(_, [:species, :nestingtype])
    DataFrames.transform(_, [:species, :nestingtype] .=> (x -> Int64.(denserank(x))) => (x -> string("num_", x)))
    groupby(_, :nestingtype)
    DataFrames.transform(_, :species => denserank => :within_nesting)
end

# Create model and prediction inputs for presence

# The length of these two vectors is used to scale up the others
num_nesting_unknown = df_species_unknown.num_nestingtype
num_region_unknown = Int64.(denserank(envs_unknown.region))

str_nesting_unknown = map(x -> x == 1 ? "burrow" : x == 2 ? "ground" : "vegetation", num_nesting_unknown)
long_nesting_unknown_num = [fill.(num_nesting_unknown, length(num_region_unknown))...;]
long_nesting_unknown_str = [fill.(str_nesting_unknown, length(num_region_unknown))...;]

# Atolls. ...
num_atoll_unknown = Int64.(denserank(envs_unknown.atoll))
long_atoll_unknown_num = [fill.(num_atoll_unknown, length(num_nesting_unknown))...;]
long_atoll_unknown_str = [fill.(envs_unknown.atoll, length(num_nesting_unknown))...;]


atoll = (
    known=(num=Int64.(denserank(envs_known.atoll)), str=envs_known.atoll),
    unknown=(num=long_atoll_unknown_num, str=long_atoll_unknown_str),
)


long_region_unknown_num = [fill.(num_region_unknown, length(num_nesting_unknown))...;]
long_region_unknown_str = [fill.(envs_unknown.region, length(num_nesting_unknown))...;]

region = (
    known=(num=Int64.(denserank(envs_known.region)), str=envs_known.region),
    unknown=(num=long_region_unknown_num, str=long_region_unknown_str),
)

num_species_unknown = [fill(df_species_unknown.num_species, length(num_region_unknown))...;]
str_species_unknown = [fill(df_species_unknown.species, length(num_region_unknown))...;]

species = (
    known=(num=Int64.(denserank(envs_known.species)), str=envs_known.species),
    unknown=(num=num_species_unknown, str=str_species_unknown)
)

unique_nesting = @chain envs_known begin
    unique(_, [:nestingtype, :species])
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    getproperty(_, :nestingtype)
    denserank(_)
    map(Int64, _)
end

nesting = (
    known=(num=Int64.(denserank(envs_known.nestingtype)), str=envs_known.nestingtype),
    unknown=(num=long_nesting_unknown_num, str=long_nesting_unknown_str),
    levels=unique_nesting,
)

presence = Float64.(envs_known.presence)

PC_known = Matrix(envs_known[:, FEATURES])
PC_unknown = Matrix(envs_unknown[!, begin:6])
long_PC_unknown = reduce(vcat, [permutedims(hcat(fill(s, length(num_nesting_unknown))...)) for s in eachslice(PC_unknown, dims=1)])

PC = (known=PC_known, unknown=long_PC_unknown)


unique_species_in_nesting = @chain envs_known begin
    unique(_, [:nestingtype, :species])
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
    map(Int64, _)
end

num_species_in_nesting = @chain envs_known begin
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    select(_, :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
end

num_species_in_nesting_unknown = df_species_unknown.within_nesting
long_num_species_in_nesting_unknown = [fill(num_species_in_nesting_unknown, length(num_region_unknown))...;]

species_in_nesting = (
    known=num_species_in_nesting,
    unknown=long_num_species_in_nesting_unknown,
    levels=unique_species_in_nesting,
)
    
end
