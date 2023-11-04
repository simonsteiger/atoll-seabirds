module PresenceVariables

export num_atoll,
       num_atoll_unknown,
       str_atoll_unknown,
       num_species,
       num_species_unknown,
       presence,
       PC,
       PC_unknown,
       count_species_by_nesting,
       num_species_within_nesting,
       num_nesting,
       num_nesting_unknown,
       num_region,
       num_region_unknown,
       unique_nesting,
       unique_species_within_nesting,
       num_species_within_nesting_unknown

include("preprocess.jl")
using .Preprocess
using DataFrames, Chain, StatsBase

const FEATURES = [:PC1, :PC2, :PC3, :PC4, :PC5, :PC6]

count_species_by_nesting = @chain begin
    unique(envs_known, :species)
    groupby(_, :nestingtype)
    combine(_, nrow)
    getproperty(_, :nrow)
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
num_atoll = Int64.(denserank(envs_known.atoll))
num_atoll_unknown = Int64.(denserank(envs_unknown.atoll))
str_atoll_unknown = envs_unknown.atoll

num_region = Int64.(denserank(envs_known.region))
num_region_unknown = Int64.(denserank(envs_unknown.region))

num_species = Int64.(denserank(envs_known.species))
num_species_unknown = df_species_unknown.num_species
str_species_unknown = df_species_unknown.species

presence = Float64.(envs_known.presence)
PC = Matrix(envs_known[:, FEATURES])
PC_unknown = Matrix(envs_unknown[!, begin:6])

num_nesting = Int64.(denserank(envs_known.nestingtype))
num_nesting_unknown = df_species_unknown.num_nestingtype

unique_nesting = @chain envs_known begin
    unique(_, [:nestingtype, :species])
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    getproperty(_, :nestingtype)
    denserank(_)
    map(Int64, _)
end

unique_species_within_nesting = @chain envs_known begin
    unique(_, [:nestingtype, :species])
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
    map(Int64, _)
end

num_species_within_nesting = @chain envs_known begin
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    select(_, :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
end

num_species_within_nesting_unknown = df_species_unknown.within_nesting
    
end
