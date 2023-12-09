module CountVariables

export num_atoll_known,
       num_atoll_unknown,
       num_region_known,
       num_region_unknown,
       num_species_known,
       str_species_known,
       num_species_unknown,
       str_species_unknown,
       nbirds,
       ppres,
       PC_known,
       PC_unknown,
       num_nesting_known,
       num_nesting_unknown,
       unique_nesting_known,
       unique_nesting_unknown,
       unique_species_within_nesting_known,
       unique_species_within_nesting_unknown,
       num_species_within_nesting_known,
       num_species_within_nesting_unknown,
       count_species_by_nesting
   
include("global.jl")
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
num_atoll_known = Int64.(denserank(pop_known.atoll))
num_atoll_unknown = Int64.(denserank(pop_unknown.atoll))

num_region_known = Int64.(denserank(pop_known.region))
num_region_unknown = Int64.(denserank(pop_unknown.region))

num_species_known = Int64.(denserank(pop_known.species))
str_species_known = pop_known.species
num_species_unknown = Int64.(denserank(pop_unknown.species))
str_species_unknown = pop_unknown.species

maximum(num_species_unknown), maximum(num_species_known)

nbirds = Float64.(pop_known.nbirds)
ppres = Float64.(pop_unknown.ppres)
PC_known = Matrix(pop_known[!, FEATURES])
PC_unknown = Matrix(pop_unknown[!, FEATURES])

num_nesting_known = Int64.(denserank(pop_known.nestingtype))
num_nesting_unknown = Int64.(denserank(pop_unknown.nestingtype))

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
unique_nesting_unknown = unique_pop_var(pop_unknown, [:nestingtype, :species], :nestingtype)

unique_species_within_nesting_known = unique_pop_var(pop_known, [:nestingtype, :species], :ne_sp)
unique_species_within_nesting_unknown = unique_pop_var(pop_unknown, [:nestingtype, :species], :ne_sp)

num_species_within_nesting_known = @chain pop_known begin
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    select(_, :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
end

num_species_within_nesting_unknown = @chain pop_unknown begin
    DataFrames.transform(_, [:nestingtype, :species] => ByRow((x,y) -> "$x$y") => :ne_sp)
    select(_, :ne_sp)
    getproperty(_, :ne_sp)
    denserank(_)
end

end