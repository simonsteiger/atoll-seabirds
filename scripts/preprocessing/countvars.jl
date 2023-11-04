module CountVariables

export num_atoll,
       num_region,
       num_species,
       nbirds,
       PC,
       num_nesting,
       unique_nesting,
       unique_species_within_nesting,
       num_species_within_nesting,
       count_species_by_nesting
   
include("preprocess.jl")
using .Preprocess

using DataFrames, Chain, StatsBase

const FEATURES = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]

# Repeat for population / count model

count_species_by_nesting = @chain begin
    unique(pop_known, :species)
    groupby(_, :nestingtype)
    combine(_, nrow)
    getproperty(_, :nrow)
end

# Atolls. ...
num_atoll = Int64.(denserank(pop_known.atoll))
# num_atoll_unknown = Int64.(denserank(envs_unknown.atoll))
# str_atoll_unknown = envs_unknown.atoll

num_region = Int64.(denserank(pop_known.region))
# num_region_unknown = Int64.(denserank(envs_unknown.region))

num_species = Int64.(denserank(pop_known.species))
# num_species_unknown = df_species_unknown.num_species
# str_species_unknown = df_species_unknown.species

nbirds = Float64.(pop_known.nbirds)
PC = Matrix(pop_known[:, FEATURES])
# PC_unknown = Matrix(envs_unknown[!, begin:6])

num_nesting = Int64.(denserank(pop_known.nestingtype))
# num_nesting_unknown = df_species_unknown.num_nestingtype

unique_nesting = @chain pop_known begin
    unique(_, [:nestingtype, :species])
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    getproperty(_, :nestingtype)
    denserank(_)
    map(Int64, _)
end

unique_species_within_nesting = @chain pop_known begin
    unique(_, [:nestingtype, :species])
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
    map(Int64, _)
end

num_species_within_nesting = @chain pop_known begin
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    select(_, :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
end

# num_species_within_nesting_unknown = df_species_unknown.within_nesting

end