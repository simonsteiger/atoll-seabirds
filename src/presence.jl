module PresenceVariables

export num_atoll,
       num_atoll_unknown,
       str_atoll_unknown,
       num_species,
       str_species,
       num_species_unknown,
       presence,
       PC,
       PC_unknown,
       count_species_by_nesting,
       num_species_in_nesting,
       num_nesting,
       num_nesting_unknown,
       num_region,
       num_region_unknown,
       unique_nesting,
       unique_species_in_nesting,
       num_species_in_nesting_unknown

include("global.jl")
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

# Atolls. ...
atoll = (
    known=(num=Int64.(denserank(envs_known.atoll)), str=envs_known.atoll),
    unknown=(num=Int64.(denserank(envs_unknown.atoll)), str=envs_unknown.atoll),
)

region = (
    known=(num=Int64.(denserank(envs_known.region)), str=envs_known.region),
    unknown=(num=Int64.(denserank(envs_unknown.region)), str=envs_unknown.region),
)

species = (
    known=(num=Int64.(denserank(envs_known.species)), str=envs_known.species),
    unknown=(num=Int64.(denserank(envs_unknown.species)), str=envs_unknown.species),
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
    unknown=(num=Int64.(denserank(envs_unknown.nestingtype)), str=envs_unknown.nestingtype),
    levels=unique_nesting,
)

presence = Float64.(envs_known.presence)

PC_known = Matrix(envs_known[:, FEATURES])
PC_unknown = Matrix(envs_unknown[!, begin:6])
PC = (known=PC, unknown=PC_unknown)


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

species_in_nesting = (
    known=num_species_in_nesting,
    unknown=num_species_in_nesting_unknown,
    levels=unique_species_in_nesting,
)
    
end
